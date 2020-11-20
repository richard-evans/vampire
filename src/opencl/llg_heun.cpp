//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <vector>
#include <sstream>

#include "atoms.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "kernels.hpp"
#include "llg_heun.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL_TIME_KERNELS
#include <chrono>

// Macro to time enqueued kernels
// will leave redundant syncs so only use for testing
#define TIME(func, var, q)                                      \
   {                                                            \
      auto start = std::chrono::high_resolution_clock::now();   \
      func;                                                     \
      q.finish();                                               \
      auto end = std::chrono::high_resolution_clock::now();     \
      std::chrono::duration<double> diff = end - start;         \
      var += diff.count();                                      \
   }
#else
#define TIME(func, var, q) func
#endif

#ifdef OPENCL
namespace vcl = ::vopencl::internal;
#endif

namespace vopencl
{
   void llg_heun(void)
   {
#ifdef OPENCL
      vcl::llg::step();
#endif // OPENCL
   }
}

#ifdef OPENCL
namespace vopencl
{
   namespace internal
   {
      cl::Kernel update_nexch_spin_fields;
      cl::Kernel update_ext;

      namespace rng
      {
         cl::Kernel grng;
      }

      namespace exchange
      {
         cl::Kernel calculate_exchange;
      }

      namespace llg
      {
         cl::Kernel predictor_step;
         cl::Kernel corrector_step;
      }

      static inline void update_spin_fields(void)
      {
         TIME(vcl::kernel_call(vcl::exchange::calculate_exchange), vcl::time::mat_mul, vcl::queue);

         vcl::queue.finish();

         TIME(vcl::kernel_call(vcl::update_nexch_spin_fields), vcl::time::spin_fields, vcl::queue);
      }

      namespace llg
      {
         void step(void) noexcept
         {
            // update random numbers for external field
            // this is independent of the rest, so if possible do it on a
            // different device or just another queue
#ifdef ENABLE_MULTIPLE_DEVICES
            TIME(vcl::kernel_call(rng::grng, vcl::queue_other, vcl::global_other), vcl::time::rng, vcl::queue_other);
#else
            const cl::CommandQueue rand_q(vcl::context, vcl::default_device);
            TIME(vcl::kernel_call(rng::grng, rand_q), vcl::time::rng, rand_q);
#endif // ENABLE_MULTIPLE_DEVICES

            // make a copy of current spins for both Heun steps
            vcl::atoms::spin_array.copy_to_dev(vcl::queue, vcl::llg::spin_buffer_array);

            // update fields
            vcl::update_spin_fields();


            // update system applied field and temperature
            vcl::update_ext.setArg(5, vcl::real_t(sim::H_vec[0] * sim::H_applied));
            vcl::update_ext.setArg(6, vcl::real_t(sim::H_vec[1] * sim::H_applied));
            vcl::update_ext.setArg(7, vcl::real_t(sim::H_vec[2] * sim::H_applied));
            vcl::update_ext.setArg(8, vcl::real_t(sim::temperature));

#ifdef ENABLE_MULTIPLE_DEVICES
            // copy random numbers to main device
            if (::gpu::platform_other != ::gpu::platform) // context to context copy
            {
               // TODO: find a way to do this device to device
               std::vector<vcl::real_t> grands(vcl::rng::n_rands);

               vcl::queue_other.finish();

               vcl::queue_other.enqueueReadBuffer(::vcl::rng::grands, CL_TRUE, 0, vcl::rng::n_rands * sizeof (vcl::real_t), grands.data());
               vcl::queue_other.enqueueWriteBuffer(::vcl::rng::grands_copy, CL_TRUE, 0, vcl::rng::n_rands * sizeof (vcl::real_t), grands.data());
            }
            else if (::gpu::device_other != ::gpu::device) // same context, different device
            {
               vcl::queue_other.finish();

               vcl::queue.enqueueCopyBuffer(::vcl::rng::grands, vcl::rng::grands_copy, 0, 0, vcl::rng::n_rands * sizeof (vcl::real_t));
            }
#else
            rand_q.finish();
#endif // ENABLE_MULTIPLE_DEVICES

            vcl::queue.finish();

            // update external and dipole fields
            TIME(vcl::kernel_call(update_ext), vcl::time::external_fields, vcl::queue);
            vcl::queue.finish();
            update_dipolar_fields();
            vcl::queue.finish();

            // Heun predictor step
            TIME(vcl::kernel_call(predictor_step), vcl::time::predictor_step, vcl::queue);
            vcl::queue.finish();

            // update spin fields, external fixed
            vcl::update_spin_fields();
            vcl::queue.finish();

            // Heun corrector step
            TIME(vcl::kernel_call(corrector_step), vcl::time::corrector_step, vcl::queue);
            vcl::queue.finish();
         }
      }
   }
}
#endif // OPENCL
