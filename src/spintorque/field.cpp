//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <iostream>

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
       
      if(st::internal::enabled==false) return;


      // update magnetisations
      st::internal::update_cell_magnetisation(x_spin_array, y_spin_array, z_spin_array, atom_type_array, mu_s_array);

      // calculate spin_accumulation
      st::internal::calculate_spin_accumulation();

      // determine spin torque field and copy to atomic internal field arrays
      for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
            const int cell3 = 3*st::internal::atom_st_index[atom];
            const double i_mu_s = 1.0/(mu_s_array[atom_type_array[atom]]);
            st::internal::x_field_array[atom] = st::internal::spin_torque[cell3+0]*i_mu_s;
            st::internal::y_field_array[atom] = st::internal::spin_torque[cell3+1]*i_mu_s;
            st::internal::z_field_array[atom] = st::internal::spin_torque[cell3+2]*i_mu_s;
            /*std::cout << atom << "\t" << st::internal::x_field_array[atom] << "\t";
            std::cout << st::internal::y_field_array[atom] << "\t" << st::internal::z_field_array[atom] << "\t";
            std::cout << i_mu_s << "\t" << mu_s_array[atom_type_array[atom]] << std::endl;
            std::cin.get();*/
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
     
      if(st::internal::enabled==false) return;


      // Add spin torque fields
      for(int i=start_index; i<end_index; ++i) x_total_external_field_array[i] += st::internal::x_field_array[i];
      for(int i=start_index; i<end_index; ++i) y_total_external_field_array[i] += st::internal::y_field_array[i];
      for(int i=start_index; i<end_index; ++i) z_total_external_field_array[i] += st::internal::z_field_array[i];

      return;
   }

} // end of st namespace
