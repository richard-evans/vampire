//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "spintorque.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{

   //-----------------------------------------------------------------------------
   // Function for updating spin torque fields
   //-----------------------------------------------------------------------------
   void update_spin_torque_fields(const std::vector<double>& x_spin_array,
                                  const std::vector<double>& y_spin_array,
                                  const std::vector<double>& z_spin_array,
                                  const std::vector<int>& atom_type_array,
                                  const std::vector<double>& mu_s_array){

      // update magnetisations
      st::internal::update_cell_magnetisation(x_spin_array, y_spin_array, z_spin_array, atom_type_array, mu_s_array);

      // calculate spin_accumulation
      st::internal::calculate_spin_accumulation();

      // copy spin torque field to atomic internal field arrays
      for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
            const int cell = st::internal::atom_st_index[atom];
            st::internal::x_field_array[atom] = st::internal::spin_torque[cell+0];
            st::internal::y_field_array[atom] = st::internal::spin_torque[cell+1];
            st::internal::z_field_array[atom] = st::internal::spin_torque[cell+2];
      }

   }

   //-----------------------------------------------------------------------------
   // Function for adding atomic spin torque fields to external field array
   //-----------------------------------------------------------------------------
   void get_spin_torque_fields(std::vector<double>& x_total_external_field_array,
                               std::vector<double>& y_total_external_field_array,
                               std::vector<double>& z_total_external_field_array,
                               const int start_index,
                               const int end_index){

      // Add spin torque fields
      for(int i=start_index; i<end_index; ++i) x_total_external_field_array[i] += st::internal::x_field_array[i];
      for(int i=start_index; i<end_index; ++i) y_total_external_field_array[i] += st::internal::y_field_array[i];
      for(int i=start_index; i<end_index; ++i) z_total_external_field_array[i] += st::internal::z_field_array[i];

      return;
   }

} // end of st namespace
