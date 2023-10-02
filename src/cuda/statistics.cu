//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "internal.hpp"
#include "statistics.hpp"
#include "typedefs.hpp"

#ifdef CUDA
namespace cu = vcuda::internal;
#endif

namespace vcuda{

#ifdef CUDA

namespace stats{
   void update (){

      // If enabled use CPU to calculate statistics by copying data from GPU
      if(vcuda::internal::stats::use_cpu){

         // copy spin data to CPU
         /*
         thrust::copy(internal::atoms::x_spin_array.begin(),internal::atoms::x_spin_array.end(),::atoms::x_spin_array.begin());
         thrust::copy(internal::atoms::y_spin_array.begin(),internal::atoms::y_spin_array.end(),::atoms::y_spin_array.begin());
         thrust::copy(internal::atoms::z_spin_array.begin(),internal::atoms::z_spin_array.end(),::atoms::z_spin_array.begin());
         */

         cudaMemcpy(internal::h_x_spin_transfer_buffer, internal::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);
         cudaMemcpy(internal::h_y_spin_transfer_buffer, internal::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);
         cudaMemcpy(internal::h_z_spin_transfer_buffer, internal::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);

         std::copy(internal::h_x_spin_transfer_buffer, internal::h_x_spin_transfer_buffer + ::atoms::num_atoms, ::atoms::x_spin_array.begin());
         std::copy(internal::h_y_spin_transfer_buffer, internal::h_y_spin_transfer_buffer + ::atoms::num_atoms, ::atoms::y_spin_array.begin());
         std::copy(internal::h_z_spin_transfer_buffer, internal::h_z_spin_transfer_buffer + ::atoms::num_atoms, ::atoms::z_spin_array.begin());

         //---------------------------------------------------------------
         // call cpu statistics functions
         //---------------------------------------------------------------

         // update magnetization statistics
         if(::stats::calculate_system_magnetization)          ::stats::system_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_material_magnetization)        ::stats::material_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_height_magnetization)          ::stats::height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_material_height_magnetization) ::stats::material_height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_grain_magnetization)           ::stats::grain_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_material_grain_magnetization)  ::stats::material_grain_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
         if(::stats::calculate_material_grain_height_magnetization) ::stats::material_grain_height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);

         // update torque statistics
         //if(stats::calculate_system_torque)          stats::system_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);
         //if(stats::calculate_grain_torque)           stats::grain_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);
         //if(stats::calculate_material_torque)        stats::material_torque.calculate_torque(sx,sy,sz,bxs,bys,bzs,bxe,bye,bze,mm);

         // update specific heat statistics
         //if(stats::calculate_system_specific_heat)         stats::system_specific_heat.calculate(stats::system_energy.get_total_energy());
         //if(stats::calculate_grain_specific_heat)          stats::grain_specific_heat.calculate(stats::grain_energy.get_total_energy());
         //if(stats::calculate_material_specific_heat)       stats::material_specific_heat.calculate(stats::material_energy.get_total_energy());

         // standard deviation in time-step
         if(::stats::calculate_material_standard_deviation)  stats::material_standard_deviation.update(stats::system_magnetization.get_magnetization());

         // update susceptibility statistics
         if(::stats::calculate_system_susceptibility)        stats::system_susceptibility.calculate(stats::system_magnetization.get_magnetization());
         if(::stats::calculate_grain_susceptibility)         stats::grain_susceptibility.calculate(stats::grain_magnetization.get_magnetization());
         if(::stats::calculate_material_susceptibility)      stats::material_susceptibility.calculate(stats::material_magnetization.get_magnetization());

         // return before doing the GPU version
         return;
      }


      // increase the counter
      cu::stats::counter++;

   }

   void get (){

      // If CPU stats calculation do nothing
      if(vcuda::internal::stats::use_cpu) return;
   }

   void reset (){

      // reset magnetization statistics
      if(vcuda::internal::stats::use_cpu){

         // reset energy statistics
         //if(stats::calculate_system_energy)                 stats::system_energy.reset_averages();
         //if(stats::calculate_grain_energy)                  stats::grain_energy.reset_averages();
         //if(stats::calculate_material_energy)               stats::material_energy.reset_averages();

         // reset magnetization statistics
         if(::stats::calculate_system_magnetization)          stats::system_magnetization.reset_magnetization_averages();
         if(::stats::calculate_grain_magnetization)           stats::grain_magnetization.reset_magnetization_averages();
         if(::stats::calculate_material_magnetization)        stats::material_magnetization.reset_magnetization_averages();
         if(::stats::calculate_material_grain_magnetization)  stats::material_grain_magnetization.reset_magnetization_averages();
         if(::stats::calculate_height_magnetization)          stats::height_magnetization.reset_magnetization_averages();
         if(::stats::calculate_material_height_magnetization) stats::material_height_magnetization.reset_magnetization_averages();
         if(::stats::calculate_material_grain_height_magnetization) stats::material_grain_height_magnetization.reset_magnetization_averages();

         // update torque statistics
         //if(stats::calculate_system_torque)          stats::system_torque.reset_torque_averages();
         //if(stats::calculate_grain_torque)           stats::grain_torque.reset_torque_averages();
         //if(stats::calculate_material_torque)        stats::material_torque.reset_torque_averages();

         // standard deviation in time-step
         if(::stats::calculate_material_standard_deviation)     stats::material_standard_deviation.reset_averages();

         // reset specific_heat statistics
         //if(stats::calculate_system_specific_heat)   stats::system_specific_heat.reset_averages();
         //if(stats::calculate_grain_specific_heat)    stats::grain_specific_heat.reset_averages();
         //if(stats::calculate_material_specific_heat) stats::material_specific_heat.reset_averages();

         // reset susceptibility statistics
         if(::stats::calculate_system_susceptibility)   stats::system_susceptibility.reset_averages();
         if(::stats::calculate_grain_susceptibility)    stats::grain_susceptibility.reset_averages();
         if(::stats::calculate_material_susceptibility) stats::material_susceptibility.reset_averages();

         // reset binder cumulant statistics
         //if(stats::calculate_system_binder_cumulant)   stats::system_binder_cumulant.reset_averages();
         //if(stats::calculate_material_binder_cumulant) stats::material_binder_cumulant.reset_averages();

         return;
      }
   }

} /* stats */

   namespace internal
   {
      namespace stats
      {
      } /* stats */
   } /* internal */

#endif

} // end of namespace cuda
