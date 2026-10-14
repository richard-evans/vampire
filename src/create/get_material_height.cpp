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

   //-------------------------------------------------------------------------------
   // Function to extract whether a material fills voided space in a
   // substructure (create::internal::mp_t::sub_fill)
   //-------------------------------------------------------------------------------
   bool get_material_sub_fill(const int material){
   	return internal::mp[material].sub_fill;
   }

   //-------------------------------------------------------------------------------
   // Function to extract a material's associated unit cell id
   // (create::internal::mp_t::unit_cell_category)
   //-------------------------------------------------------------------------------
   int get_material_unit_cell_category(const int material){
   	return internal::mp[material].unit_cell_category;
   }

   //-------------------------------------------------------------------------------
   // Function to extract a material's nucleation height for the grain
   // substructure path (create::internal::mp_t::voronoi_grain_substructure_nucleation_height)
   //-------------------------------------------------------------------------------
   double get_material_substructure_nucleation_height(const int material){
   	return internal::mp[material].voronoi_grain_substructure_nucleation_height;
   }

} // end of namespace create
