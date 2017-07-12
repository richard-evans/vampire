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
   // Function to calculate energy from anisotropy tensors
   //
   // Note - tensors are defined with units of Tesla and a double calcellation
   //        of minus (-) signs, one from the Hamiltonian and the other from the
   //        derivative in the field calculation.
   //
   //---------------------------------------------------------------------------
   double single_spin_energy(const int atom, const int imaterial, const double sx, const double sy, const double sz, const double temperature){

      // get index for tensor
      const unsigned int index = 9*atom;

      // variable to add energies
      double energy = 0.0;

      // second order tensor
      if(internal::enable_second_order_tensor){
         energy = energy - (  sx * internal::second_order_tensor[index + 0] * sx
                            + sx * internal::second_order_tensor[index + 1] * sy
                            + sx * internal::second_order_tensor[index + 2] * sz
                            + sy * internal::second_order_tensor[index + 3] * sx
                            + sy * internal::second_order_tensor[index + 4] * sy
                            + sy * internal::second_order_tensor[index + 5] * sz
                            + sz * internal::second_order_tensor[index + 6] * sx
                            + sz * internal::second_order_tensor[index + 7] * sy
                            + sz * internal::second_order_tensor[index + 8] * sz);
      }

      // fourth order tensor
      if(internal::enable_fourth_order_tensor){

         // unroll constants
         const double sx2 = sx*sx;
         const double sy2 = sy*sy;
         const double sz2 = sz*sz;

         energy = energy - (  sx2 * internal::fourth_order_tensor[index + 0] * sx2
                            + sx2 * internal::fourth_order_tensor[index + 1] * sy2
                            + sx2 * internal::fourth_order_tensor[index + 2] * sz2
                            + sy2 * internal::fourth_order_tensor[index + 3] * sx2
                            + sy2 * internal::fourth_order_tensor[index + 4] * sy2
                            + sy2 * internal::fourth_order_tensor[index + 5] * sz2
                            + sz2 * internal::fourth_order_tensor[index + 6] * sx2
                            + sz2 * internal::fourth_order_tensor[index + 7] * sy2
                            + sz2 * internal::fourth_order_tensor[index + 8] * sz2);
      }

      // optionally calclulate lattice anisotropy fields
      if(internal::enable_lattice_anisotropy){
         energy += internal::spin_lattice_anisotropy_energy(imaterial, sx, sy, sz, temperature);
      }

      return energy;

   }

} // end of anisotropy namespace
