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
#include "spintorque.hpp"
#include "vmpi.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Funtion to calculate the magnetisation of all cells
      //-----------------------------------------------------------------------------
      void update_cell_magnetisation(const std::vector<double>& x_spin_array,
                                     const std::vector<double>& y_spin_array,
                                     const std::vector<double>& z_spin_array,
                                     const std::vector<int>& atom_type_array,
                                     const std::vector<double>& mu_s_array){

         // calculate number of cells and number of elements
         const std::vector<int>::size_type num_elements = st::internal::m.size();

         // reset cell magnetisations to zero
         for(int i=0; i<num_elements; ++i) st::internal::m[i]=0.0;

         // calulate total moment in each cell
         for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
            const int cell = st::internal::atom_st_index[atom];
            const int material = atom_type_array[atom];
            const double mus = mu_s_array[material];

            st::internal::m[3*cell+0] += x_spin_array[atom]*mus;
            st::internal::m[3*cell+1] += y_spin_array[atom]*mus;
            st::internal::m[3*cell+2] += z_spin_array[atom]*mus;

            //calculate average mus of each ST cell
            const double natom = st::internal::cell_natom[cell];
            const double i_natom = 1.0/natom;
            st::internal::cell_mus[cell] += mus*i_natom;



         }

         #ifdef MPICF
            // Add all microcell magnetisations on all nodes
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::m[0],st::internal::m.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
      }

   } // end of namespace internal
} // end of namespace st
