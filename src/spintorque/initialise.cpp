//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <complex>
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
   //if(!st::internal::enabled) return;
    if(st::internal::enabled==false) return;

   // output informative message
   zlog << zTs() << "Initialising data structures for spin torque calculation." << std::endl;

   //-------------------------------------------------------------------------------------
   // Determine transformation between x,y,z in actual and spin torque coordinate system
   //-------------------------------------------------------------------------------------
   int stx=0; // indices for x,y,z in the spin torque coordinate system (default z)
   int sty=1;
   int stz=2;
   if(st::internal::current_direction==0){
      stx=2; // c[stx] = c[2] = atom_x // current direction
      sty=1; // c[sty] = c[1] = atom_y
      stz=0; // c[stz] = c[0] = atom_z
   }
   else if(st::internal::current_direction==1){
      stx=0;// c[stx] = c[0] = atom_x
      sty=2;// c[sty] = c[2] = atom_y // current direction
      stz=1;// c[stz] = c[1] = atom_z
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
   st::internal::diffusion.resize(array_size); /// diffusion constant Do
   st::internal::sd_exchange.resize(array_size); /// diffusion constant Do
   st::internal::a.resize(array_size); // a parameter for spin accumulation
   st::internal::b.resize(array_size); // b parameter for spin accumulation
   st::internal::coeff_ast.resize(array_size);
   st::internal::coeff_nast.resize(array_size);
   st::internal::cell_natom.resize(array_size);


   const int three_vec_array_size = 3*array_size;
   st::internal::pos.resize(three_vec_array_size); /// microcell position
   st::internal::m.resize(three_vec_array_size); // magnetisation
   st::internal::j.resize(three_vec_array_size); // spin current
   st::internal::sa.resize(three_vec_array_size); // spin accumulation
   st::internal::spin_torque.resize(three_vec_array_size); // spin torque
   st::internal::ast.resize(three_vec_array_size); // adiabatic spin torque
   st::internal::nast.resize(three_vec_array_size); // non-adiabatic spin torque
   st::internal::total_ST.resize(three_vec_array_size); // non-adiabatic spin torque



   //---------------------------------------------------
   // Noi Initialise j,sa, st, ast, nast here?
   //---------------------------------------------------
   for(int cell = 0; cell < array_size; ++cell){
      st::internal::sa[3*cell+0] = 0.0;
      st::internal::sa[3*cell+1] = 0.0;
      st::internal::sa[3*cell+2] = 1.0;
      st::internal::j [3*cell+0] = 0.0;
      st::internal::j [3*cell+1] = 0.0;
      st::internal::j [3*cell+2] = 0.0;
   }

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

   // Allocate space for 3D supercell array (ST coordinate system)
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

   // define array to store atom-microcell associations
   st::internal::atom_st_index.resize(num_local_atoms);

   // Determine number of cells in x,y,z (ST coordinate system)
   const int d[3]={ncx,ncy,ncz};
   const double cs[3] = {st::internal::micro_cell_size, st::internal::micro_cell_size, st::internal::micro_cell_thickness}; // cell size

   // Assign atoms to cells
   for(int atom=0;atom<num_local_atoms;atom++){
      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[stx]=atom_coords_x[atom]+0.0001;
      c[sty]=atom_coords_y[atom]+0.0001;
      c[stz]=atom_coords_z[atom]+0.0001;
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
            std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
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

   // optionally output base microcell data
   st::internal::output_base_microcell_data();

   return;
}

namespace internal{

   //--------------------------------------------------------------------------------
   // Function to determine spin torque properties from atomic material properties
   //--------------------------------------------------------------------------------
   void set_microcell_properties(const std::vector<int>& atom_type_array, const int num_local_atoms){
      //-------------------------------------------------------
      // Determine microcell properties from atomic properties
      //-------------------------------------------------------
      st::internal::default_properties.beta_cond = 0.11;
      st::internal::default_properties.beta_diff = 0.36;
      st::internal::default_properties.sa_infinity = 1.0e8;
      st::internal::default_properties.lambda_sdl = 100.0e-9; // m
      st::internal::default_properties.diffusion = 0.0001; //Angstroms^2/s
      st::internal::default_properties.sd_exchange = 1.6e-21; //Joule


      // Temporary array to hold number of atoms in each cell for averaging
      std::vector<double> count(st::internal::beta_cond.size(),0.0);

      // loop over all atoms
      for(int atom=0;atom<num_local_atoms;atom++){

         // get material type
         int mat = atom_type_array[atom];

         // get microcell id
         int id = st::internal::atom_st_index[atom];

         // determine atomic properties
         double beta_cond = st::internal::mp.at(mat).beta_cond; // beta
         double beta_diff = st::internal::mp.at(mat).beta_diff; // beta_prime
         double sa_infinity = st::internal::mp.at(mat).sa_infinity;
         double lambda_sdl = st::internal::mp.at(mat).lambda_sdl;
         double diffusion = st::internal::mp.at(mat).diffusion;
         double sd_exchange = st::internal::mp.at(mat).sd_exchange;

         //add atomic properties to microcells
         st::internal::beta_cond.at(id) += beta_cond;
         st::internal::beta_diff.at(id) += beta_diff;
         st::internal::sa_infinity.at(id) += sa_infinity;
         st::internal::lambda_sdl.at(id) += lambda_sdl;
         st::internal::diffusion.at(id) += diffusion;
         st::internal::sd_exchange.at(id) += sd_exchange;
         count.at(id) += 1.0;
      }

      // reduce microcell properties on all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_cond[0],   st::internal::beta_cond.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_diff[0],   st::internal::beta_diff.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sa_infinity[0], st::internal::sa_infinity.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::lambda_sdl[0],  st::internal::lambda_sdl.size(),  MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::diffusion[0],   st::internal::diffusion.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sd_exchange[0], st::internal::sd_exchange.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &count[0],                     count.size(),                     MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      // Calculate average (mean) spin torque parameters
      for(size_t cell=0; cell<beta_cond.size(); ++cell){
         const double nat = count.at(cell);

          st::internal::cell_natom[cell] = nat;



         // check for zero atoms in cell
         if(nat>0.0001){
            st::internal::beta_cond.at(cell)   /= nat;
            st::internal::beta_diff.at(cell)   /= nat;
            st::internal::sa_infinity.at(cell) /= nat;
            st::internal::lambda_sdl.at(cell)  /= nat;
            st::internal::diffusion.at(cell)   /= nat;
            st::internal::sd_exchange.at(cell) /= nat;
         }

         else{
            st::internal::beta_cond.at(cell)   = st::internal::default_properties.beta_cond;
            st::internal::beta_diff.at(cell)   = st::internal::default_properties.beta_diff;
            st::internal::sa_infinity.at(cell) = st::internal::default_properties.sa_infinity;
            st::internal::lambda_sdl.at(cell)  = st::internal::default_properties.lambda_sdl;
            st::internal::diffusion.at(cell)   = st::internal::default_properties.diffusion;
            st::internal::sd_exchange.at(cell) = st::internal::default_properties.sd_exchange;
         }
      }



      // Determine a and b parameters
      const double hbar = 1.05457162e-34;
      for(size_t cell=0; cell<beta_cond.size(); ++cell){

         const double B  = st::internal::beta_cond[cell];
         const double Bp = st::internal::beta_diff[cell];
         const double lambda_sdl = st::internal::lambda_sdl[cell];
         const double Do = st::internal::diffusion[cell];
         const double Jsd = st::internal::sd_exchange[cell];

         const double BBp = 1.0/sqrt(1.0-B*Bp);
         const double lambda_sf = lambda_sdl*BBp;
         const double lambda_j = sqrt(2.0*hbar*Do/Jsd); // Angstroms
         const double lambda_sf2 = lambda_sf*lambda_sf;
         const double lambda_j2 = lambda_j*lambda_j;

         std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
         std::complex<double> inv_lplus = sqrt(inside);

         st::internal::a[cell] =  real(inv_lplus);
         st::internal::b[cell] = -imag(inv_lplus);

      }

      return;
   }
}

} // end of st namespace
