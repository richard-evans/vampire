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

   namespace stats
   {
         void update ()
         {

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

               // call cpu statistics functions
               if(::stats::calculate_system_magnetization)          ::stats::system_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_material_magnetization)        ::stats::material_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_height_magnetization)          ::stats::height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_material_height_magnetization) ::stats::material_height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);

               // return before doing the GPU version
               return;
            }


            // increase the counter
            cu::stats::counter++;

         }

         void get ()
         {

            // If CPU stats calculation do nothing
            if(vcuda::internal::stats::use_cpu) return;
         }

         void reset ()
         {
             // reset magnetization statistics
             if(vcuda::internal::stats::use_cpu){
                 if(::stats::calculate_system_magnetization)          ::stats::system_magnetization.reset_magnetization_averages();
                 if(::stats::calculate_material_magnetization)        ::stats::material_magnetization.reset_magnetization_averages();
                 if(::stats::calculate_height_magnetization)          ::stats::height_magnetization.reset_magnetization_averages();
                 if(::stats::calculate_material_height_magnetization) ::stats::material_height_magnetization.reset_magnetization_averages();
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
