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

   int num_atoms; // Number of atoms for statistic purposes

   bool calculate_system_energy                 = false;
   bool calculate_grain_energy                  = false;
   bool calculate_material_energy               = false;
   
   bool calculate_system_sld_energy                 = false;
   bool calculate_grain_sld_energy                  = false;
   bool calculate_material_sld_energy               = false;

   bool calculate_system_magnetization          = true;
   bool calculate_grain_magnetization           = false;
   bool calculate_material_magnetization        = false;
   bool calculate_material_grain_magnetization  = false;
   bool calculate_height_magnetization          = false;
   bool calculate_material_height_magnetization = false;

   bool calculate_system_torque                 = false;
   bool calculate_grain_torque                  = false;
   bool calculate_material_torque               = false;

   bool calculate_system_spin_temp                 = false;
   bool calculate_grain_spin_temp                  = false;
   bool calculate_material_spin_temp              = false;
   
   bool calculate_system_lattice_temp                 = false;
   bool calculate_grain_lattice_temp                  = false;
   bool calculate_material_lattice_temp              = false;
 
   bool calculate_system_specific_heat          = false;
   bool calculate_grain_specific_heat           = false;
   bool calculate_material_specific_heat        = false;

   bool calculate_material_standard_deviation   = false;

   bool calculate_system_susceptibility         = false;
   bool calculate_grain_susceptibility          = false;
   bool calculate_material_susceptibility       = false;

   energy_statistic_t system_energy("s");
   energy_statistic_t grain_energy("g");
   energy_statistic_t material_energy("m");
   
   sld_energy_statistic_t system_sld_energy("s");
   sld_energy_statistic_t grain_sld_energy("g");
   sld_energy_statistic_t material_sld_energy("m");

   magnetization_statistic_t system_magnetization("s");
   magnetization_statistic_t grain_magnetization("g");
   magnetization_statistic_t material_magnetization("m");
   magnetization_statistic_t material_grain_magnetization("mg");
   magnetization_statistic_t height_magnetization("h");
   magnetization_statistic_t material_height_magnetization("mh");

   torque_statistic_t system_torque("s");
   torque_statistic_t grain_torque("g");
   torque_statistic_t material_torque("m");
   
   spin_temp_statistic_t system_spin_temp("s");
   spin_temp_statistic_t grain_spin_temp("g");
   spin_temp_statistic_t material_spin_temp("m");
   
   lattice_temp_statistic_t system_lattice_temp("s");
   lattice_temp_statistic_t grain_lattice_temp("g");
   lattice_temp_statistic_t material_lattice_temp("m");

   specific_heat_statistic_t system_specific_heat("s");
   specific_heat_statistic_t grain_specific_heat("s");
   specific_heat_statistic_t material_specific_heat("m");

   standard_deviation_statistic_t material_standard_deviation("m");

   susceptibility_statistic_t system_susceptibility("s");
   susceptibility_statistic_t grain_susceptibility("g");
   susceptibility_statistic_t material_susceptibility("m");

   //-----------------------------------------------------------------------------
   // Shared variables used for statistics calculation
   //-----------------------------------------------------------------------------
   namespace internal{

   } // end of internal namespace
} // end of stats namespace
