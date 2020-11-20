//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo Meo 2020. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "create.hpp"
#include "internal.hpp"

namespace create{

   //-------------------------------------------------------------------------------
   // Functions to extract min fractional height of material
   //-------------------------------------------------------------------------------
   double get_material_height_min(const int material){
   	return internal::mp[material].min;
   }

   //-------------------------------------------------------------------------------
   // Functions to extract max fractional height of material
   //-------------------------------------------------------------------------------
   double get_material_height_max(const int material){
   	return internal::mp[material].max;
   }

} // end of namespace create
