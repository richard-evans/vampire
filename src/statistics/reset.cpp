//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "stats.hpp"

namespace stats{

   //------------------------------------------------------------------------------------------------------
   // Function to reset required statistics classes
   //------------------------------------------------------------------------------------------------------
   void reset(){

      // Check for GPU acceleration and reset statistics on device
      if(gpu::acceleration){
         gpu::stats::reset();
      }
      else{
         // reset energy statistics
         if(stats::calculate_system_energy)                 stats::system_energy.reset_averages();
         if(stats::calculate_grain_energy)                  stats::grain_energy.reset_averages();
         if(stats::calculate_material_energy)               stats::material_energy.reset_averages();

         // reset magnetization statistics
         if(stats::calculate_system_magnetization)          stats::system_magnetization.reset_magnetization_averages();
         if(stats::calculate_grain_magnetization)           stats::grain_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_magnetization)        stats::material_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_grain_magnetization)  stats::material_grain_magnetization.reset_magnetization_averages();
         if(stats::calculate_height_magnetization)          stats::height_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_height_magnetization) stats::material_height_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_grain_height_magnetization) stats::material_grain_height_magnetization.reset_magnetization_averages();

         // update torque statistics
         if(stats::calculate_system_torque)          stats::system_torque.reset_torque_averages();
         if(stats::calculate_grain_torque)           stats::grain_torque.reset_torque_averages();
         if(stats::calculate_material_torque)        stats::material_torque.reset_torque_averages();

         // standard deviation in time-step
         if(stats::calculate_material_standard_deviation)     stats::material_standard_deviation.reset_averages();

         // reset specific_heat statistics
         if(stats::calculate_system_specific_heat)   stats::system_specific_heat.reset_averages();
         if(stats::calculate_grain_specific_heat)    stats::grain_specific_heat.reset_averages();
         if(stats::calculate_material_specific_heat) stats::material_specific_heat.reset_averages();

         // reset susceptibility statistics
         if(stats::calculate_system_susceptibility)   stats::system_susceptibility.reset_averages();
         if(stats::calculate_grain_susceptibility)    stats::grain_susceptibility.reset_averages();
         if(stats::calculate_material_susceptibility) stats::material_susceptibility.reset_averages();

         // reset binder cumulant statistics
         if(stats::calculate_system_binder_cumulant)   stats::system_binder_cumulant.reset_averages();
         if(stats::calculate_material_binder_cumulant) stats::material_binder_cumulant.reset_averages();

      }

      return;

   }

}
