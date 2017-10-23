//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"
#include "atoms.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy {

   double second_order_tensor_anisotropy(){

      double energy = 0;

      for (int atom=0; atom<internal::num_atoms; ++atom){

         /* get atom spin */
         double sx = atoms::x_spin_array.at(atom);
         double sy = atoms::y_spin_array.at(atom);
         double sz = atoms::z_spin_array.at(atom);

         /* first matrix product [(1x3).(3x3) => (1x3)] */
         double ux =
         sx * internal::second_order_tensor.at(atom).at(0).at(0) +
         sy * internal::second_order_tensor.at(atom).at(0).at(1) +
         sy * internal::second_order_tensor.at(atom).at(0).at(2);

         double uy =
         sx * internal::second_order_tensor.at(atom).at(1).at(0) +
         sy * internal::second_order_tensor.at(atom).at(1).at(1) +
         sy * internal::second_order_tensor.at(atom).at(1).at(2);

         double uz =
         sx * internal::second_order_tensor.at(atom).at(2).at(0) +
         sy * internal::second_order_tensor.at(atom).at(2).at(1) +
         sy * internal::second_order_tensor.at(atom).at(2).at(2);

         /* second matrix product [(1x3).(3x1) => scalar] */
         energy = ux * sx +
         uy * sy +
         uz * sz;

      }

      return energy;
   }

} // end of anisotropy namespace
