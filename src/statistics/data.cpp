//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "stats.hpp"

namespace stats{

   bool calculate_system_magnetization          = true;
   bool calculate_material_magnetization        = false;
   bool calculate_height_magnetization          = false;
   bool calculate_material_height_magnetization = false;
   bool calculate_system_susceptibility         = false;

   magnetization_statistic_t system_magnetization;
   magnetization_statistic_t material_magnetization;
   magnetization_statistic_t height_magnetization;
   magnetization_statistic_t material_height_magnetization;

   susceptibility_statistic_t system_susceptibility;

   //-----------------------------------------------------------------------------
   // Shared variables used for statistics calculation
   //-----------------------------------------------------------------------------
   namespace internal{

   } // end of internal namespace
} // end of stats namespace
