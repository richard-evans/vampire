//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //---------------------------------------------------------------------------
   // Function to calculate magnetic fields from anisotropy tensors
   //---------------------------------------------------------------------------
   void fields(std::vector<double>& spin_array_x,
               std::vector<double>& spin_array_y,
               std::vector<double>& spin_array_z,
               std::vector<int>&    type_array,
               std::vector<double>& field_array_x,
               std::vector<double>& field_array_y,
               std::vector<double>& field_array_z,
               const int start_index,
               const int end_index,
               const double temperature){

      // second order uniaxial anisotropy
      internal::uniaxial_second_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // fourth order uniaxial anisotropy
      internal::uniaxial_fourth_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // fourth order biaxial anisotropy (simple version)
      internal::biaxial_fourth_order_simple_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      //triaxial anisotropy variable basis
      internal::triaxial_second_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      //triaxial anisotropy fixed basis
      internal::triaxial_second_order_fields_fixed_basis(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // triaxial fourth order anisotropy
      internal::triaxial_fourth_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // triaxial fourth order anisotropy fixed basis
      internal::triaxial_fourth_order_fields_fixed_basis(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // sixth order uniaxial anisotropy
      internal::uniaxial_sixth_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // fourth order rotational anisotropy
      internal::rotational_fourth_order_fields_fixed_basis(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // fourth order cubic anisotropy
      internal::cubic_fourth_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // fourth order cubic anisotropy (rotated basis)
      internal::cubic_fourth_order_rotation_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // sixth order cubic anisotropy
      internal::cubic_sixth_order_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // Neel anisotropy
      internal::neel_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index);

      // lattice anisotropy
      internal::lattice_fields(spin_array_x, spin_array_y, spin_array_z, type_array, field_array_x, field_array_y, field_array_z, start_index, end_index, temperature);

      return;

   }

} // end of anisotropy namespace
