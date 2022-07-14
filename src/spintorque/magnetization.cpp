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
#include <iostream>
#include "material.hpp"


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

         st::internal::magx_mat.resize(mp::num_materials); // array size of magnetisation of material
         st::internal::magy_mat.resize(mp::num_materials);
         st::internal::magz_mat.resize(mp::num_materials);


         // reset cell magnetisations to zero
         for(int i=0; i<mp::num_materials; ++i) st::internal::magx_mat[i]=0.0;
         for(int i=0; i<mp::num_materials; ++i) st::internal::magy_mat[i]=0.0;
         for(int i=0; i<mp::num_materials; ++i) st::internal::magz_mat[i]=0.0;
         for(unsigned int i=0; i<num_elements; ++i) st::internal::m[i]=0.0;

         // calculate total moment in each cell
         for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
            const int cell = st::internal::atom_st_index[atom];
            const int material = atom_type_array[atom];
            const double mus = mu_s_array[material];

            st::internal::m[3*cell+0] += x_spin_array[atom]*mus;
            st::internal::m[3*cell+1] += y_spin_array[atom]*mus;
            st::internal::m[3*cell+2] += z_spin_array[atom]*mus;

            //calculate magnetisation of each material
            st::internal::magx_mat[material] += x_spin_array[atom];
            st::internal::magy_mat[material] += y_spin_array[atom];
            st::internal::magz_mat[material] += z_spin_array[atom];

         }

         #ifdef MPICF
            // Add all microcell magnetisations on all nodes
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::m[0],st::internal::m.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::magx_mat[0],st::internal::magx_mat.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::magy_mat[0],st::internal::magy_mat.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::magz_mat[0],st::internal::magz_mat.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         //calculate the normalised magnetisation of each material
         for(int mat=0; mat<mp::num_materials; ++mat) {
            double norm_mag = sqrt(st::internal::magx_mat[mat]*st::internal::magx_mat[mat]
                           +st::internal::magy_mat[mat]*st::internal::magy_mat[mat]
                           +st::internal::magz_mat[mat]*st::internal::magz_mat[mat]);

            if(norm_mag>1.0e-8){
		st::internal::magx_mat[mat] /= norm_mag;
		st::internal::magy_mat[mat] /= norm_mag;
		st::internal::magz_mat[mat] /= norm_mag;
            }
         }
      }

   } // end of namespace internal
} // end of namespace st
