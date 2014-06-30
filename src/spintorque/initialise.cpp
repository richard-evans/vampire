//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <vector>

// Vampire headers
#include "errors.hpp"
#include "spintorque.hpp"
#include "vio.hpp"

// Spin Torque headers
#include "internal.hpp"

//-----------------------------------------------------------------------------
// Forward function declarations
//-----------------------------------------------------------------------------
namespace st{
   namespace internal{
      void set_microcell_properties(const std::vector<int>& atom_type_array, const int num_local_atoms);
   }
}

namespace st{

//-----------------------------------------------------------------------------
// Function for initialising spin torque data structure and variables
//-----------------------------------------------------------------------------
void initialise(const double system_dimensions_x,
                const double system_dimensions_y,
                const double system_dimensions_z,
                const std::vector<double>& atom_coords_x,
                const std::vector<double>& atom_coords_y,
                const std::vector<double>& atom_coords_z,
                const std::vector<int>& atom_type_array,
                const int num_local_atoms){

   //-------------------------------------------------------------------------------------
   // Check for spin torque calculation enabled, if not do nothing
   //-------------------------------------------------------------------------------------
   if(!st::internal::enabled) return;

   // output informative message
   zlog << zTs() << "Initialising data structures for spin torque calculation." << std::endl;

   //-------------------------------------------------------------------------------------
   // Determine transformation between x,y,z in actual and spin torque coordinate system
   //-------------------------------------------------------------------------------------
   int stx=0; // indices for x,y,z in the spin torque coordinate system (default z)
   int sty=1;
   int stz=2;
   if(st::internal::current_direction==0){
      stx=1;
      sty=2;
      stz=0;
   }
   else if(st::internal::current_direction==1){
      stx=0;
      sty=2;
      stz=1;
   }

   //-------------------------------------------------------------------------------------
   // Calculate number of stacks and microcells
   //-------------------------------------------------------------------------------------

   // Make a small array for system dimensions
   double system_dimensions[3]={system_dimensions_x,system_dimensions_y,system_dimensions_z};

   // determine number of cells in each stack  (global)
   st::internal::num_microcells_per_stack = 1+ceil((system_dimensions[stz]+0.01)/st::internal::micro_cell_thickness);

   // determine number of stacks in x and y (global)
   st::internal::num_x_stacks = ceil((system_dimensions[stx]+0.01)/st::internal::micro_cell_size);
   st::internal::num_y_stacks = ceil((system_dimensions[sty]+0.01)/st::internal::micro_cell_size);

   // determine total number of stacks
   st::internal::num_stacks = st::internal::num_x_stacks*st::internal::num_y_stacks;

   // allocate array to store index of first element of stack
   st::internal::stack_index.resize(st::internal::num_stacks);
   

   //-------------------------------------------------------------------------------------
   // allocate microcell data
   //-------------------------------------------------------------------------------------
   const int array_size = st::internal::num_stacks*st::internal::num_microcells_per_stack;
   st::internal::beta_cond.resize(array_size); /// spin polarisation (conductivity)
   st::internal::beta_diff.resize(array_size); /// spin polarisation (diffusion)
   st::internal::sa_infinity.resize(array_size); /// intrinsic spin accumulation
   st::internal::lambda_sdl.resize(array_size); /// spin diffusion length
   const int three_vec_array_size = 3*array_size;
   st::internal::pos.resize(three_vec_array_size); /// microcell position
   st::internal::m.resize(three_vec_array_size); // magnetisation
   st::internal::j.resize(three_vec_array_size); // spin current
   st::internal::sa.resize(three_vec_array_size); // spin accumulation
   st::internal::spin_torque.resize(three_vec_array_size); // spin torque
   st::internal::ast.resize(three_vec_array_size); // adiabatic spin torque
   st::internal::nast.resize(three_vec_array_size); // non-adiabatic spin torque

   //---------------------------------------------------
   // Noi Initialise j,sa, st, ast, nast here?
   //---------------------------------------------------

   //---------------------------------------------------
   // Determine which atoms belong to which stacks
   //---------------------------------------------------
   {
   int ncx = st::internal::num_x_stacks; // temporary variables for readability
   int ncy = st::internal::num_y_stacks;
   int ncz = st::internal::num_microcells_per_stack;

   // Set cell and stack counters
   int cell=0;
   int stack=0;

   // Declare array for create space for 3D supercell array
   std::vector<std::vector<std::vector<int> > > supercell_array;
   supercell_array.resize(ncx);
   for(int i=0;i<ncx;++i){
      supercell_array[i].resize(ncy);
      for(int j=0;j<ncy;++j){
         supercell_array[i][j].resize(ncz);
         // set starting cell for each stack
         st::internal::stack_index[stack]=cell;
         // increment stack counter
         stack++;
         // store cell coordinates
         for(int k=0; k<ncz; ++k){
            // associate cell with position i,j,k
            supercell_array[i][j][k]=cell;
            // save ijk coordinates as microcell positions
            st::internal::pos.at(3*cell+0)=i;
            st::internal::pos.at(3*cell+1)=j;
            st::internal::pos.at(3*cell+2)=k;
            // increment cell number
            cell++;
         }
      }
   }

   // slightly offset atomic coordinates to prevent fence post problem
   double atom_offset[3]={0.01,0.01,0.01};

   // define array to store atom-microcell associations
   st::internal::atom_st_index.resize(num_local_atoms);

   // Determine number of cells in x,y,z
   const int d[3]={ncx,ncy,ncz};
   const double cs[3] = {st::internal::micro_cell_size, st::internal::micro_cell_size, st::internal::micro_cell_thickness}; // cell size

   // Assign atoms to cells                                                                                                                  
   for(int atom=0;atom<num_local_atoms;atom++){
      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[stx]=atom_coords_x[atom]+0.01;
      c[sty]=atom_coords_y[atom]+0.01;
      c[stz]=atom_coords_z[atom]+0.01;
      int scc[3]={0,0,0}; // super cell coordinates
      for(int i=0;i<3;i++){
         // Determine supercell coordinates for atom (rounding down)
         scc[i]=int(c[i]/cs[i]);
         // Always check cell in range
         if(scc[i]<0 || scc[i]>= d[i]){
            terminaltextcolor(RED);
            std::cerr << "Error - atom out of supercell range in spin torque microcell calculation!" << std::endl;
            terminaltextcolor(WHITE);
            #ifdef MPICF
            terminaltextcolor(RED);
            std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
            terminaltextcolor(WHITE);
            #endif 
            terminaltextcolor(RED);
            std::cerr << "\tAtom number:      " << atom << std::endl;
            std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
            std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
            std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      // If no error for range then assign atom to cell.
      st::internal::atom_st_index[atom]=supercell_array[scc[0]][scc[1]][scc[2]+1]; // move cells up by one in z
   }
   } // end of supercell assignment of atoms

   //-------------------------------------------------------
   // Determine microcell properties from atomic properties
   //-------------------------------------------------------
   st::internal::set_microcell_properties(atom_type_array, num_local_atoms);

   //-------------------------------------------------------
   // Save value of local num atoms and resize field arrays
   //-------------------------------------------------------
   st::internal::num_local_atoms = num_local_atoms;
   st::internal::x_field_array.resize(num_local_atoms); // arrays to store atomic spin torque field
   st::internal::y_field_array.resize(num_local_atoms);
   st::internal::z_field_array.resize(num_local_atoms);

   st::internal::output_microcell_data();
   return; 
}

namespace internal{

   //--------------------------------------------------------------------------------
   // Function to determine spin torque properties from atomic material properties
   //--------------------------------------------------------------------------------
   void set_microcell_properties(const std::vector<int>& atom_type_array, const int num_local_atoms){

      st::internal::mp.resize(1);
      st::internal::mp.at(0).beta_cond=1.0;
      st::internal::mp.at(0).beta_diff=2.0;
      st::internal::mp.at(0).sa_infinity=3.0;
      st::internal::mp.at(0).lambda_sdl=4.0;

      //-------------------------------------------------------
      // Determine microcell properties from atomic properties
      //-------------------------------------------------------
      st::internal::default_properties.beta_cond = 5.0;
      st::internal::default_properties.beta_diff = 5.0;
      st::internal::default_properties.sa_infinity = 5.0;
      st::internal::default_properties.lambda_sdl = 5.0;

      // Temporary array to hold number of atoms in each cell for averaging
      std::vector<double> count(st::internal::beta_cond.size(),0.0);

      // loop over all atoms
      for(int atom=0;atom<num_local_atoms;atom++){

         // get material type
         int mat = atom_type_array[atom];

         // get microcell id
         int id = st::internal::atom_st_index[atom];

         // determine atomic properties
         double beta_cond = st::internal::mp.at(mat).beta_cond;
         double beta_diff = st::internal::mp.at(mat).beta_diff;
         double sa_infinity = st::internal::mp.at(mat).sa_infinity;
         double lambda_sdl = st::internal::mp.at(mat).lambda_sdl;

         //add atomic properties to microcells
         st::internal::beta_cond.at(id) += beta_cond;
         st::internal::beta_diff.at(id) += beta_diff;
         st::internal::sa_infinity.at(id) += sa_infinity;
         st::internal::lambda_sdl.at(id) += lambda_sdl;
         count.at(id) += 1.0;
      }

      // reduce microcell properties on all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_cond[0],st::internal::beta_cond.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_diff[0],st::internal::beta_diff.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sa_infinity[0],st::internal::sa_infinity.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::lambda_sdl[0],st::internal::lambda_sdl.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &count[0], count.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      // Calculate average (mean) spin torque parameters
      for(int cell=0; cell<beta_cond.size(); ++cell){
         const double nat = count.at(cell);
         // check for zero atoms in cell
         if(nat>0.0001){
            st::internal::beta_cond.at(cell)/=nat;
            st::internal::beta_diff.at(cell)/=nat;
            st::internal::sa_infinity.at(cell)/=nat;
            st::internal::lambda_sdl.at(cell)/=nat;
         }
         else{
            st::internal::beta_cond.at(cell)=st::internal::default_properties.beta_cond;
            st::internal::beta_diff.at(cell)=st::internal::default_properties.beta_cond;
            st::internal::sa_infinity.at(cell)=st::internal::default_properties.sa_infinity;
            st::internal::lambda_sdl.at(cell)=st::internal::default_properties.lambda_sdl;
         }
      }

      return;
   }
}

} // end of st namespace

