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
#include "atoms.hpp"
#include "gpu.hpp"
#include "sim.hpp"
#include "stats.hpp"

namespace stats{

   //-----------------------------------------------------------------------------
   // Shared variables used for statistics calculation
   //-----------------------------------------------------------------------------
   namespace internal{

      //------------------------------------------------------------------------------------------------------
      // Function to update required statistics classes
      //------------------------------------------------------------------------------------------------------
      void update(const std::vector<double>& sx, // spin unit vector
                  const std::vector<double>& sy,
                  const std::vector<double>& sz,
                  const std::vector<double>& bxs, // spin fields
                  const std::vector<double>& bys,
                  const std::vector<double>& bzs,
                  const std::vector<double>& bxe, // external fields
                  const std::vector<double>& bye,
                  const std::vector<double>& bze,
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
            if(stats::calculate_grain_energy)                  stats::grain_energy.calculate(sx, sy, sz, mm, mat, temperature);
            if(stats::calculate_material_energy)               stats::material_energy.calculate(sx, sy, sz, mm, mat, temperature);

            // update magnetization statistics
            if(stats::calculate_system_magnetization)          stats::system_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_grain_magnetization)           stats::grain_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_material_magnetization)        stats::material_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_material_grain_magnetization)  stats::material_grain_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_height_magnetization)          stats::height_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_material_height_magnetization) stats::material_height_magnetization.calculate_magnetization(sx,sy,sz,mm);
            if(stats::calculate_material_grain_height_magnetization) stats::material_grain_height_magnetization.calculate_magnetization(sx,sy,sz,mm);

            // update torque statistics
            if(stats::calculate_system_torque)          stats::system_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);
            if(stats::calculate_grain_torque)           stats::grain_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);
            if(stats::calculate_material_torque)        stats::material_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);

            // update specific heat statistics
            if(stats::calculate_system_specific_heat)         stats::system_specific_heat.calculate(stats::system_energy.get_total_energy());
            if(stats::calculate_grain_specific_heat)          stats::grain_specific_heat.calculate(stats::grain_energy.get_total_energy());
            if(stats::calculate_material_specific_heat)       stats::material_specific_heat.calculate(stats::material_energy.get_total_energy());

            // standard deviation in time-step
            if(stats::calculate_material_standard_deviation)  stats::material_standard_deviation.update(stats::system_magnetization.get_magnetization());

            // update susceptibility statistics
            if(stats::calculate_system_susceptibility)        stats::system_susceptibility.calculate(stats::system_magnetization.get_magnetization());
            if(stats::calculate_grain_susceptibility)         stats::grain_susceptibility.calculate(stats::grain_magnetization.get_magnetization());
            if(stats::calculate_material_susceptibility)      stats::material_susceptibility.calculate(stats::material_magnetization.get_magnetization());

            // update binder cumulant statistics
            if(stats::calculate_system_binder_cumulant)         stats::system_binder_cumulant.calculate(stats::system_magnetization.get_magnetization());
            if(stats::calculate_material_binder_cumulant)       stats::material_binder_cumulant.calculate(stats::material_magnetization.get_magnetization());

            // update spin length statistics
            if(stats::calculate_system_spin_length)           stats::system_spin_length.calculate_spin_length(sx,sy,sz);
            if(stats::calculate_material_spin_length)         stats::material_spin_length.calculate_spin_length(sx,sy,sz);
            if(stats::calculate_height_spin_length)           stats::height_spin_length.calculate_spin_length(sx,sy,sz);

         }

         return;

      }

   } // end of internal namespace

   //------------------------------------------------------------------------------------------------------
   // Wrapper function to update required statistics classes
   //------------------------------------------------------------------------------------------------------
   void update(){

      // call actual function, picking up arguments directly from namespace header files
      stats::internal::update(atoms::x_spin_array, 				  		atoms::y_spin_array, 				    atoms::z_spin_array,
   					            atoms::x_total_spin_field_array,     atoms::y_total_spin_field_array, 	 atoms::z_total_spin_field_array,
   					            atoms::x_total_external_field_array, atoms::y_total_external_field_array, atoms::z_total_external_field_array,
   					            atoms::m_spin_array, 					   atoms::type_array, 						 sim::temperature);

      return;

   }


} // end of stats namespace
