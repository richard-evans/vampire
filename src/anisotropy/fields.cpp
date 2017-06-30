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
   //
   // Hx = -dE/dSx = - 2      (  Txx Sx  +  Txy Sy  +  Txz Sz  )
   //                - 4 Sx   ( Txx Sx^2 + Txy Sy^2 + Txz Sz^2 )
   //                - 6 Sx^2 ( Txx Sx^3 + Txy Sy^3 + Txz Sz^3 )
   //
   // Hy = -dE/dSy = - 2      (  Tyx Sx  +  Tyy Sy  +  Tyz Sz  )
   //                - 4 Sy   ( Tyx Sx^2 + Tyy Sy^2 + Tyz Sz^2 )
   //                - 6 Sy^2 ( Tyx Sx^3 + Tyy Sy^3 + Tyz Sz^3 )
   //
   // Hz = -dE/dSz = - 2      (  Tzx Sx  +  Tzy Sy  +  Tzz Sz  )
   //                - 4 Sz   ( Tzx Sx^2 + Tzy Sy^2 + Tzz Sz^2 )
   //                - 6 Sz^2 ( Tzx Sx^3 + Tzy Sy^3 + Tzz Sz^3 )
   //
   // Note - tensors are defined with units of Tesla and a double calcellation
   //        of minus (-) signs, one from the Hamiltonian and the other from the
   //        derivative in the field calculation.
   //
   //---------------------------------------------------------------------------
   void fields(std::vector<double>& spin_array_x,
               std::vector<double>& spin_array_y,
               std::vector<double>& spin_array_z,
               std::vector<double>& field_array_x,
               std::vector<double>& field_array_y,
               std::vector<double>& field_array_z,
               const int start_index,
               const int end_index){

      // Loop over all atoms between start and end index
      for(int atom = start_index; atom<end_index; atom++){

         const double sx = spin_array_x[atom]; // store spin direction in temporary variables
         const double sy = spin_array_y[atom];
         const double sz = spin_array_z[atom];

         const unsigned int index = 9*atom; // get atom index in tensor array

         double hx = 2.0 * ( internal::second_order_tensor[index + 0] * sx +
                             internal::second_order_tensor[index + 1] * sy +
                             internal::second_order_tensor[index + 2] * sz );

         double hy = 2.0 * ( internal::second_order_tensor[index + 3] * sx +
                             internal::second_order_tensor[index + 4] * sy +
                             internal::second_order_tensor[index + 5] * sz );

         double hz = 2.0 * ( internal::second_order_tensor[index + 6] * sx +
                             internal::second_order_tensor[index + 7] * sy +
                             internal::second_order_tensor[index + 8] * sz );

         // store net field
         field_array_x[atom] += hx;
         field_array_y[atom] += hy;
         field_array_z[atom] += hz;

      }

      return;

   }

} // end of anisotropy namespace
