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

#define TIME(func, var)                                         \
   {                                                            \
      auto start = std::chrono::high_resolution_clock::now();   \
      func;                                                     \
      vcl::queue.finish();                                      \
      auto end = std::chrono::high_resolution_clock::now();     \
      std::chrono::duration<double> diff = end - start;         \
      var += diff.count();                                      \
   }
#else
#define TIME(func, var) func
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

      static inline void update_spin_fields(void)
      {
         //vcl::total_spin_field_array.zero_buffer();

         TIME(vcl::kernel_call(vcl::update_nexch_spin_fields), vcl::time::spin_fields);

         vcl::queue.finish();

         TIME(vcl::kernel_call(vcl::exchange::calculate_exchange), vcl::time::mat_mul);
      }

      namespace llg
      {
         cl::Kernel predictor_step;
         cl::Kernel corrector_step;

         void step(void) noexcept
         {
            // make a copy of current spins for both Heun steps
            vcl::atoms::spin_array.copy_to_dev(vcl::queue, vcl::llg::spin_buffer_array);

            // update fields
            vcl::update_spin_fields();

            // update random numbers for external field
            TIME(vcl::kernel_call(rng::grng), vcl::time::rng);
            vcl::queue.finish();

            // update system applied field and temperature
            vcl::update_ext.setArg(5, vcl::real_t(sim::H_vec[0] * sim::H_applied));
            vcl::update_ext.setArg(6, vcl::real_t(sim::H_vec[1] * sim::H_applied));
            vcl::update_ext.setArg(7, vcl::real_t(sim::H_vec[2] * sim::H_applied));
            vcl::update_ext.setArg(8, vcl::real_t(sim::temperature));

            // update external and dipole fields
            TIME(vcl::kernel_call(update_ext), vcl::time::external_fields);
            vcl::queue.finish();
            update_dipolar_fields();
            vcl::queue.finish();

            // Heun predictor step
            TIME(vcl::kernel_call(predictor_step), vcl::time::predictor_step);
            vcl::queue.finish();

            // update spin fields, external fixed
            vcl::update_spin_fields();
            vcl::queue.finish();

            // Heun corrector step
            TIME(vcl::kernel_call(corrector_step), vcl::time::corrector_step);
            vcl::queue.finish();
         }
      }
   }
}
#endif // OPENCL
