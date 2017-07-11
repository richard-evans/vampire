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

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //--------------------------------------------------------------------------------
   // Function to get second order uniaxial anisotropy constant for a given material
   //--------------------------------------------------------------------------------
   double get_anisotropy_constant(const int material){
      return internal::mp[material].ku2;
   }

} // end of anisotropy namespace
