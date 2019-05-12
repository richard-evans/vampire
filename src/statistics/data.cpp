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

   bool calculate_system_energy                 = false;
   bool calculate_material_energy               = false;
   bool calculate_system_magnetization          = true;
   bool calculate_material_magnetization        = false;
   bool calculate_height_magnetization          = false;
   bool calculate_material_height_magnetization = false;
   bool calculate_system_specific_heat          = false;
   bool calculate_material_specific_heat        = false;
   bool calculate_material_standard_deviation   = false;
   bool calculate_system_susceptibility         = false;
   bool calculate_material_susceptibility       = false;

   energy_statistic_t system_energy;
   energy_statistic_t material_energy;

   //torque_statistic_t system_torque;

   magnetization_statistic_t system_magnetization;
   magnetization_statistic_t material_magnetization;
   magnetization_statistic_t height_magnetization;
   magnetization_statistic_t material_height_magnetization;

   specific_heat_statistic_t system_specific_heat;
   specific_heat_statistic_t material_specific_heat;

   
   standard_deviation_statistic_t material_standard_deviation;
   susceptibility_statistic_t system_susceptibility;
   susceptibility_statistic_t material_susceptibility;

   //-----------------------------------------------------------------------------
   // Shared variables used for statistics calculation
   //-----------------------------------------------------------------------------
   namespace internal{

   } // end of internal namespace
} // end of stats namespace
