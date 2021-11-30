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
   // Function to calculate anisotropy energy
   //---------------------------------------------------------------------------
   double single_spin_energy(const int atom, const int mat, const double sx, const double sy, const double sz, const double temperature){

      // variable to add energies
      double energy = 0.0;

      // if not enabled then do nothing
      if(internal::enable_uniaxial_second_order)         energy += internal::uniaxial_second_order_energy(atom, mat, sx, sy, sz);
      if(internal::enable_uniaxial_fourth_order)         energy += internal::uniaxial_fourth_order_energy(atom, mat, sx, sy, sz);
      if(internal::enable_biaxial_fourth_order_simple)   energy += internal::biaxial_fourth_order_simple_energy(atom, mat, sx, sy, sz);
      if(internal::enable_uniaxial_sixth_order)          energy += internal::uniaxial_sixth_order_energy (atom, mat, sx, sy, sz);

      if(internal::enable_cubic_fourth_order)            energy += internal::cubic_fourth_order_energy(atom, mat, sx, sy, sz);
      if(internal::enable_cubic_fourth_order_rotation)   energy += internal::cubic_fourth_order_rotation_energy(atom, mat, sx, sy, sz);
      if(internal::enable_cubic_sixth_order)             energy += internal::cubic_sixth_order_energy (atom, mat, sx, sy, sz);

      if(internal::enable_fourth_order_rotational)       energy += internal::rotational_fourth_order_energy_fixed_basis(atom, mat, sx, sy, sz);

      if(internal::enable_neel_anisotropy)               energy += internal::neel_energy(atom, mat, sx, sy, sz);
      if(internal::enable_lattice_anisotropy)            energy += internal::lattice_energy(atom, mat, sx, sy, sz, temperature);

      if(internal::enable_triaxial_anisotropy)           energy += internal::triaxial_second_order_energy_fixed_basis(atom, mat, sx, sy, sz);
      if(internal::enable_triaxial_anisotropy_rotated)   energy += internal::triaxial_second_order_energy(atom, mat, sx, sy, sz);
      if(internal::enable_triaxial_fourth_order)         energy += internal::triaxial_fourth_order_energy_fixed_basis(atom, mat, sx, sy, sz);
      if(internal::enable_triaxial_fourth_order_rotated) energy += internal::triaxial_fourth_order_energy(atom, mat, sx, sy, sz);

      return energy;

   }

} // end of anisotropy namespace
