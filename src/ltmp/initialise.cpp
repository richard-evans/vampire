//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "ltmp.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vio.hpp"

// Local temperature pulse headers
#include "internal.hpp"

struct uvec{
  int i;
  int j;
  int k;
};

namespace ltmp{

//---------------------------------------------------------------------------------
// Function for initialising local temperature pulse data structures and variables
//---------------------------------------------------------------------------------
void initialise(const double system_dimensions_x,
                const double system_dimensions_y,
                const double system_dimensions_z,
                const std::vector<double>& atom_coords_x,
                const std::vector<double>& atom_coords_y,
                const std::vector<double>& atom_coords_z,
                const std::vector<int>& atom_type_array,
                const int num_local_atoms,
                const double starting_temperature,
                const double pump_power,
                const double pump_time,
                const double TTG,
                const double TTCe,
                const double TTCl,
                const double dt,
                const double Tmin,
                const double Tmax
               ){

   //-------------------------------------------------------------------------------------
   // Check for local temperature pulse calculation enabled, if not do nothing
   //-------------------------------------------------------------------------------------
   if(!ltmp::internal::enabled) return;

   // output informative message
   zlog << zTs() << "Initialising data structures for local temperature calculation." << std::endl;

   // check for prior initialisation
   if(ltmp::internal::initialised){
      zlog << zTs() << "Warning: Localised temperature calculation already initialised. Continuing." << std::endl;
      return;
   }

   //-------------------------------------------------------------------------------------
   // Set simulation constants
   //-------------------------------------------------------------------------------------
   ltmp::internal::pump_power=pump_power; // laser pump power
   ltmp::internal::pump_time=pump_time; // laser pump time
   ltmp::internal::TTG=TTG;  // electron-lattice coupling constant
   ltmp::internal::TTCe=TTCe; // electron heat capacity (T=0)
   ltmp::internal::TTCl=TTCl; // lattice heat capcity
   ltmp::internal::dt=dt; // time step
   ltmp::internal::minimum_temperature = Tmin; // minimum temperature for temperature gradient
   ltmp::internal::maximum_temperature = Tmax; // maximum temperature for temperature gradient

   //-------------------------------------------------------------------------------------
   // Calculate number of microcells
   //-------------------------------------------------------------------------------------

   // determine number of stacks in x and y (global)
   int dx =  ceil((system_dimensions_x+0.01)/ltmp::internal::micro_cell_size);
   int dy =  ceil((system_dimensions_y+0.01)/ltmp::internal::micro_cell_size);
   int dz =  ceil((system_dimensions_z+0.01)/ltmp::internal::micro_cell_size);

   // determine total number of stacks
   if(     ltmp::internal::lateral_discretisation == true  && ltmp::internal::vertical_discretisation == true ){
      ltmp::internal::num_cells = dx*dy*dz;
   }
   else if(ltmp::internal::lateral_discretisation == true  && ltmp::internal::vertical_discretisation == false){
      ltmp::internal::num_cells = dx*dy;
      dz = 1;
   }
   else if(ltmp::internal::lateral_discretisation == false && ltmp::internal::vertical_discretisation == true ){
      ltmp::internal::num_cells = dz;
      dx = 1;
      dy = 1;
   }
   // no discretisation - disable local thermal fields and return
   else{
      zlog << zTs() << "No discretisation (lateral/vertical) specified for local temperature calculation. Disabling local temperature calculation." << std::endl;
      ltmp::internal::enabled = false;
      return;
   }

   //-------------------------------------------------------------------------------------
   // Allocate microcell data and initialise starting temperature (Teq)
   //-------------------------------------------------------------------------------------
   const double sqrt_starting_temperature = sqrt(starting_temperature);
   ltmp::internal::root_temperature_array.resize(2*ltmp::internal::num_cells,sqrt_starting_temperature);
   ltmp::internal::cell_position_array.resize(3*ltmp::internal::num_cells);

   //---------------------------------------------------
   // Determine which atoms belong to which cell
   //---------------------------------------------------
   {
   int ncx = dx; // temporary variables for readability
   int ncy = dy;
   int ncz = dz;

   // Set cell and stack counters
   int cell=0;

   // allocate temporary array for neighbour list calculation
   std::vector<uvec> cell_list;
   cell_list.reserve(dx*dy*dz);

   // Allocate space for 3D supercell array (ST coordinate system)
   std::vector<std::vector<std::vector<int> > > supercell_array;
   supercell_array.resize(ncx);
   for(int i=0;i<ncx;++i){
      supercell_array[i].resize(ncy);
      for(int j=0;j<ncy;++j){
         supercell_array[i][j].resize(ncz);
         // store cell coordinates
         for(int k=0; k<ncz; ++k){

            // associate cell with position i,j,k
            supercell_array[i][j][k]=cell;

            // reverse association for neighbour calculation (store i,j,k)
            uvec tmp;
            tmp.i = i;
            tmp.j = j;
            tmp.k = k;
            cell_list.push_back(tmp);

            // save ijk coordinates as microcell positions
            ltmp::internal::cell_position_array.at(3*cell+0)=double(i)*ltmp::internal::micro_cell_size;
            ltmp::internal::cell_position_array.at(3*cell+1)=double(j)*ltmp::internal::micro_cell_size;
            ltmp::internal::cell_position_array.at(3*cell+2)=double(k)*ltmp::internal::micro_cell_size;

            // increment cell number
            cell++;
         }
      }
   }

   // define array to store atom-microcell associations
   ltmp::internal::atom_temperature_index.resize(num_local_atoms);

   // Determine number of cells in x,y,z (ST coordinate system)
   const int d[3]={ncx,ncy,ncz};
   const double cs[3] = {ltmp::internal::micro_cell_size, ltmp::internal::micro_cell_size, ltmp::internal::micro_cell_size}; // cell size

   // Assign atoms to cells
   for(int atom=0;atom<num_local_atoms;atom++){
      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[0]=atom_coords_x[atom]+0.0001;
      c[1]=atom_coords_y[atom]+0.0001;
      c[2]=atom_coords_z[atom]+0.0001;
      int scc[3]={0,0,0}; // super cell coordinates
      // Determine supercell coordinates for atom (rounding down)
      if(     ltmp::internal::lateral_discretisation == true  && ltmp::internal::vertical_discretisation == true ){
         scc[0]=int(c[0]/cs[0]);
         scc[1]=int(c[1]/cs[1]);
         scc[2]=int(c[2]/cs[2]);
      }
      else if(     ltmp::internal::lateral_discretisation == true  && ltmp::internal::vertical_discretisation == false ){
         scc[0]=int(c[0]/cs[0]);
         scc[1]=int(c[1]/cs[1]);
         scc[2]=0;
      }
      else if(     ltmp::internal::lateral_discretisation == false  && ltmp::internal::vertical_discretisation == true ){
         scc[0]=0;
         scc[1]=0;
         scc[2]=int(c[2]/cs[2]);
      }
      for(int i=0;i<3;i++){
         // Always check cell in range
         if(scc[i]<0 || scc[i]>= d[i]){
            terminaltextcolor(RED);
            std::cerr << "Error - atom out of supercell range in local temperature microcell calculation!" << std::endl;
            terminaltextcolor(WHITE);
            #ifdef MPICF
            terminaltextcolor(RED);
            std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
            terminaltextcolor(WHITE);
            #endif
            terminaltextcolor(RED);
            std::cerr << "\tAtom number:      " << atom << std::endl;
            std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
            std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
            std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
            std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      // If no error for range then assign atom to cell
      int cell = supercell_array[scc[0]][scc[1]][scc[2]];
      // Now determine whether atom couples to electron or phonon temperature
      int mat = atom_type_array[atom];
      if(mp::material[mat].couple_to_phonon_temperature){
         ltmp::internal::atom_temperature_index[atom] = 2*cell+1;
      }
      else{
         ltmp::internal::atom_temperature_index[atom]= 2*cell+0;
      }
   }

   //-------------------------------------------------------------------------------------
   // Determine number of microcells computed locally
   //-------------------------------------------------------------------------------------
   int num_local_cells = ltmp::internal::num_cells; // parallelise!!

   ltmp::internal::delta_temperature_array.resize(2*num_local_cells);
   ltmp::internal::attenuation_array.resize(num_local_cells);

   // calculate laser position
   double laser_x = system_dimensions_x*0.5;
   double laser_y = system_dimensions_y*0.5;

   // Only parallelise for lateral discretisation
   //if(ltmp::internal::lateral_discretisation){
     // int num_cells_per_cpu = ltmp::internal::num_cells/vmpi::num_processors;
     // ltmp::internal::my_first_cell=vmpi::my_rank
   //}
   //MPI_Gather(void* send_data, int send_count, MPI_Datatype send_datatype,
     //      void* recv_data, int recv_count, MPI_Datatype recv_datatype,
     //      int root, MPI_Comm communicator)

   // determine if profile is taken from file
   const bool profile_file=ltmp::absorption_profile.is_set();

   // calculate interpolation for absorption profile
   ltmp::absorption_profile.set_interpolation_table();

   //--------------------------------------------------------------------------------------
   // calculate attenuation for each cell depending on lateral and vertical discretisation
   //--------------------------------------------------------------------------------------
   for(int cell = 0; cell < num_local_cells; ++cell){
      double rx = ltmp::internal::cell_position_array[3*cell+0]-laser_x;
      double ry = ltmp::internal::cell_position_array[3*cell+1]-laser_y;
      double r_sq = rx*rx+ry*ry;
      double z = system_dimensions_z - ltmp::internal::cell_position_array[3*cell+2];
      double pre = 4.0*log(2.0);
      double vattn = 1.0;
      if(ltmp::internal::vertical_discretisation){
         if(profile_file) vattn = ltmp::absorption_profile.get_absorption_constant(z);
         else vattn = exp(-z/ltmp::internal::penetration_depth); // vertical attenuation
         // Check for gradient and if so overwrite with linear profile
         if(ltmp::internal::gradient) vattn = ltmp::internal::cell_position_array[3*cell+2]/system_dimensions_z;
      }
      double lattn = 1.0;
      if(ltmp::internal::lateral_discretisation){
         lattn = exp(-r_sq*pre/(ltmp::internal::laser_spot_size*ltmp::internal::laser_spot_size));
         // Check for gradient and if so overwrite with linear profile in x
         if(ltmp::internal::gradient) lattn = ltmp::internal::cell_position_array[3*cell+0]/system_dimensions_x;
      }
      // determine attenuation
      ltmp::internal::attenuation_array[cell] = vattn*lattn;

   }

   //-------------------------------------------------------
   // Determine list of cell interactions for heat transfer
   //-------------------------------------------------------
   ltmp::internal::cell_neighbour_start_index.resize(num_local_cells);
   ltmp::internal::cell_neighbour_end_index.resize(num_local_cells);
   ltmp::internal::cell_neighbour_list.reserve(3*num_local_cells);

   int index_counter = 0;
   // loop over all cells and determine neighbouring cells
   for(int cell = 0; cell < num_local_cells; ++cell){

      const int i = cell_list[cell].i;
      const int j = cell_list[cell].j;
      const int k = cell_list[cell].k;

      // set starting index
      ltmp::internal::cell_neighbour_start_index[cell]=index_counter;

      if(i+1 >=0 && i+1 < dx){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i+1][j][k]); index_counter++;}
      if(i-1 >=0 && i-1 < dx){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i-1][j][k]); index_counter++;}
      if(j+1 >=0 && j+1 < dy){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i][j+1][k]); index_counter++;}
      if(j-1 >=0 && j-1 < dy){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i][j-1][k]); index_counter++;}
      if(k+1 >=0 && k+1 < dz){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i][j][k+1]); index_counter++;}
      if(k-1 >=0 && k-1 < dz){ ltmp::internal::cell_neighbour_list.push_back(supercell_array[i][j][k-1]); index_counter++;}

      // set end index
      ltmp::internal::cell_neighbour_end_index[cell]=index_counter;

      // Print cell neighbour list (Debugging)
      //std::cout << cell << "\t|\t" << i << "\t" << j << "\t" << k << "\t|\t";
      //for(int id=ltmp::internal::cell_neighbour_start_index[cell]; id<ltmp::internal::cell_neighbour_end_index[cell]; ++id) std::cout << ltmp::internal::cell_neighbour_list[id] << "\t";
      //std::cout << std::endl;

   }

   } // end of supercell assignment of atoms

   //-------------------------------------------------------
   // Save value of local num atoms and resize field arrays
   //-------------------------------------------------------
   ltmp::internal::num_local_atoms = num_local_atoms;
   ltmp::internal::x_field_array.resize(num_local_atoms); // arrays to store atomic spin torque field
   ltmp::internal::y_field_array.resize(num_local_atoms);
   ltmp::internal::z_field_array.resize(num_local_atoms);

   //-------------------------------------------------------
   // Unroll thermal field prefactor for all atoms
   //-------------------------------------------------------
   ltmp::internal::atom_sigma.resize(num_local_atoms);
   for(int atom=0; atom<num_local_atoms; ++atom){
      ltmp::internal::atom_sigma[atom] = mp::material[atom_type_array[atom]].H_th_sigma;
   }

   //------------------------------------------------------------------
   // Optionally unroll temperature rescaling prefactors for all atoms
   //------------------------------------------------------------------
   // Determine if rescaling is needed (slower performance) (if Tc > 0)
   for(int mat=0; mat<mp::num_materials; mat++) if(mp::material[mat].temperature_rescaling_Tc>0.0) ltmp::internal::temperature_rescaling=true;

   if(ltmp::internal::temperature_rescaling){
      ltmp::internal::atom_rescaling_root_Tc.resize(num_local_atoms);
      ltmp::internal::atom_rescaling_alpha.resize(num_local_atoms);
      for(int atom=0; atom<num_local_atoms; ++atom){
         ltmp::internal::atom_rescaling_root_Tc[atom] = sqrt(mp::material[atom_type_array[atom]].temperature_rescaling_Tc);
         ltmp::internal::atom_rescaling_alpha[atom] = mp::material[atom_type_array[atom]].temperature_rescaling_alpha;
      }
   }

   // optionally output temperature cell data
   if(ltmp::internal::output_microcell_data){
      ltmp::internal::write_microcell_data();
      // initial output file for vertical temperature profile
      if(!ltmp::internal::lateral_discretisation && ltmp::internal::vertical_discretisation) ltmp::internal::open_vertical_temperature_profile_file();
      // initial output file for lateral temperature profile
      if(ltmp::internal::lateral_discretisation && !ltmp::internal::vertical_discretisation) ltmp::internal::open_lateral_temperature_profile_file();
   }

   // Set initialised flag
   ltmp::internal::initialised = true;

   return;
}

} // end of st namespace
