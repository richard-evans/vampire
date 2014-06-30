//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Funtion to calculate the spin accumulation and spin torque
      //-----------------------------------------------------------------------------
      void calculate_spin_accumulation(){

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (st::internal::spin_torque.begin(),st::internal::spin_torque.end(),0.0);

         // loop over all 1D stacks (in parallel)
         for(int stack=0; stack <num_stacks; ++stack){
            // determine starting cell in stack
            const int idx = stack_index[stack];
            // loop over all cells in stack after first (idx+1)
            for(int cell=idx+1; cell<idx+num_microcells_per_stack; ++cell){

               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // caculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               // copy array values to temporaries for readability

               // current cell magnetisations
               const double mx = st::internal::m[cellx];
               const double my = st::internal::m[celly];
               const double mz = st::internal::m[cellz];

               // previous cell magnetisations
               const double pmx = st::internal::m[pcellx];
               const double pmy = st::internal::m[pcelly];
               const double pmz = st::internal::m[pcellz];

               //std::cout << mx << "\t" << my << "\t" << mz << std::endl;
               //... Calculate spin torque...
               // double stx = ...
               // double sty = ...
               // double stz = ...

               // maybe some functions?
               // update_j();

               // st::internal::spin_torque[cell+0] = stx;
               // st::internal::spin_torque[cell+1] = sty;
               // st::internal::spin_torque[cell+2] = stz;

            } // end of cell loop
         } // end of stack loop

         #ifdef MPICF
            // Reduce all microcell spin torques on all nodes
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_torque[0],st::internal::spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         return;
      }

   } // end of namespace internal
} // end of namespace st
