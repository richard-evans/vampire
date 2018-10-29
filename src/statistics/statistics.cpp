//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "gpu.hpp"
#include "stats.hpp"
#include "vmpi.hpp"

namespace stats{

   //------------------------------------------------------------------------------------------------------
   // Function to update required statistics classes
   //------------------------------------------------------------------------------------------------------
   void update(const std::vector<double>& sx, // spin unit vector
               const std::vector<double>& sy,
               const std::vector<double>& sz,
               const std::vector<double>& mm,
               const std::vector<int>& mat,
               const double temperature
            ){

      // Check for GPU acceleration and update statistics on device
      if(gpu::acceleration){
         gpu::stats::update();
      }
      else{
         // update energy statistics
         if(stats::calculate_system_energy)                 stats::system_energy.calculate(sx, sy, sz, mm, mat, temperature);
         if(stats::calculate_material_energy)               stats::material_energy.calculate(sx, sy, sz, mm, mat, temperature);

         // update magnetization statistics
         if(stats::calculate_system_magnetization)          stats::system_magnetization.calculate_magnetization(sx,sy,sz,mm);
         if(stats::calculate_material_magnetization)        stats::material_magnetization.calculate_magnetization(sx,sy,sz,mm);
         if(stats::calculate_height_magnetization)          stats::height_magnetization.calculate_magnetization(sx,sy,sz,mm);
         if(stats::calculate_material_height_magnetization) stats::material_height_magnetization.calculate_magnetization(sx,sy,sz,mm);

         // update specific heat statistics
         if(stats::calculate_system_specific_heat)         stats::system_specific_heat.calculate(stats::system_energy.get_total_energy());
         if(stats::calculate_material_specific_heat)       stats::material_specific_heat.calculate(stats::material_energy.get_total_energy());

         // update susceptibility statistics
         if(stats::calculate_system_susceptibility)         stats::system_susceptibility.calculate(stats::system_magnetization.get_magnetization());
         if(stats::calculate_material_susceptibility)       stats::material_susceptibility.calculate(stats::material_magnetization.get_magnetization());

      }

      return;

   }

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
         if(stats::calculate_material_energy)               stats::material_energy.reset_averages();

         // reset magnetization statistics
         if(stats::calculate_system_magnetization)          stats::system_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_magnetization)        stats::material_magnetization.reset_magnetization_averages();
         if(stats::calculate_height_magnetization)          stats::height_magnetization.reset_magnetization_averages();
         if(stats::calculate_material_height_magnetization) stats::material_height_magnetization.reset_magnetization_averages();

         // reset specific_heat statistics
         if(stats::calculate_system_specific_heat) stats::system_specific_heat.reset_averages();
         if(stats::calculate_material_specific_heat) stats::material_specific_heat.reset_averages();

         // reset susceptibility statistics
         if(stats::calculate_system_susceptibility) stats::system_susceptibility.reset_averages();
         if(stats::calculate_material_susceptibility) stats::material_susceptibility.reset_averages();

      }

      return;

   }

}
