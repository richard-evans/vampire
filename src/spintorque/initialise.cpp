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
#include "atoms.hpp"
#include "errors.hpp"
#include "vio.hpp"

namespace st{
   
//-----------------------------------------------------------------------------
// Function for initialising spin torque data structure and variables
//-----------------------------------------------------------------------------
void initialise(){
      
   // determine transformation between x,y,z in actual and spin torque coordinate system
   int stx=0; // indices for x,y,z in the spin torque coordinate system (default z)
   int sty=1;
   int stz=2;
   if(st::current_direction==0){
      stx=1;
      sty=2;
      stz=0;
   }
   else if(st::current_direction==1){
      stx=0;
      sty=2;
      stz=1;
   }
   
   // determine number of cells in each stack  (global)
   st::num_microcells_per_stack = 1+ceil((cs::system_dimensions[stz]+0.01)/st::micro_cell_thickness);
   
   // determine number of stacks in x and y (global)
   st::num_x_stacks = ceil((cs::system_dimensions[stx]+0.01)/st::micro_cell_size);
   st::num_y_stacks = ceil((cs::system_dimensions[sty]+0.01)/st::micro_cell_size);
   
   // determine total number of stacks
   st::num_stacks = num_x_stacks*num_y_stacks;

   // allocate array to store index of first element of stack
   st::stack_index.resize(num_stacks);
   
   
   
   // allocate microcell data
   const int array_size = st::num_stacks*st::num_microcells_per_stack;
   st::beta_cond.resize(array_size); /// spin polarisation (conductivity)
   st::beta_diff.resize(array_size); /// spin polarisation (diffusion)
   st::sa_infinity.resize(array_size); /// intrinsic spin accumulation
   st::lambda_sdl.resize(array_size); /// spin diffusion length
   const int three_vec_array_size = 3*array_size;
   st::pos.resize(three_vec_array_size); /// microcell position
   st::m.resize(three_vec_array_size); // magnetisation
   st::j.resize(three_vec_array_size); // spin current
   st::sa.resize(three_vec_array_size); // spin accumulation
   st::st.resize(three_vec_array_size); // spin torque
   st::ast.resize(three_vec_array_size); // adiabatic spin torque
   st::nast.resize(three_vec_array_size); // non-adiabatic spin torque      

   //---------------------------------------------------
   // Determine which atoms belong to which stacks
   //--------------------------------------------------
   {
   int ncx = st::num_x_stacks; // temporary variables for readability
   int ncy = st::num_y_stacks;
   int ncz = st::num_microcells_per_stack;

   // Set cell counter
   int cell=0;

   // Declare array for create space for 3D supercell array
   std::vector<std::vector<std::vector<int> > > supercell_array;
   supercell_array.resize(ncx);
   for(int i=0;i<ncx;++i){
      supercell_array[i].resize(ncy);
      for(int j=0;j<ncy;++j){
         supercell_array[i][j].resize(ncz);
         for(int k=0; k<ncz; ++k){
            supercell_array[i][j][k]=cell;
            cell++;
         }
      }
   }

   // slightly offset atomic coordinates to prevent fence post problem
   double atom_offset[3]={0.01,0.01,0.01};

   // For MPI version, only add local atoms
   #ifdef MPICF
      int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
   #else
      int num_local_atoms = atoms::num_atoms;
   #endif

   // define array to store atom-microcell associations
   st::atom_st_index.resize(num_local_atoms);

   // Determine number of cells in x,y,z
   const int d[3]={ncx,ncy,ncz};
   const double cs[3] = {st::micro_cell_size, st::micro_cell_size, st::micro_cell_thickness}; // cell size

   // Assign atoms to cells                                                                                                                  
   for(int atom=0;atom<num_local_atoms;atom++){
      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[stx]=atoms::x_coord_array[atom]+0.01;
      c[sty]=atoms::y_coord_array[atom]+0.01;
      c[stz]=atoms::z_coord_array[atom]+0.01;
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
      st::atom_st_index[atom]=supercell_array[scc[0]][scc[1]][scc[2]+1]; // move cells up by one in z
   }

   } // end of supercell assignment of atoms

   // reduce mean properties lambda sdl on all CPUs for all materials

   return; 
}   


      
   
   
//-----------------------------------------------------------------------------
// calculate microcell magnetisation
//-----------------------------------------------------------------------------

// loop over all atoms and add s to microcells

// all reduce to have consistent m on all CPUs
/*   for(int stack=0; stack <num_stacks; stack++){
   const int idx = start_index[stack];
   for(int cell=idx; cell<idx+num_microcells_per_stack; ++cell){
      double mx = st::mx[cell];
      //...
   }   
   
}*/

//-----------------------------------------------------------------------------
// Calculate ST field
//-----------------------------------------------------------------------------
   
// divide stacks among CPUs
   
} // end of st namespace

