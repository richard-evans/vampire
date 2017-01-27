#include <vector>
#include <sstream>

#include "atoms.hpp"
#include "material.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "kernels.hpp"
#include "llg_heun.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

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
         vcl::total_spin_field_array.zero_buffer();
         vcl::kernel_call(vcl::update_nexch_spin_fields);

         vcl::queue.finish();

         vcl::kernel_call(vcl::exchange::calculate_exchange);
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
            vcl::kernel_call(rng::grng);

            vcl::queue.finish();

            // update external and dipole fields
            vcl::kernel_call(update_ext);
            update_dipolar_fields();

            vcl::queue.finish();

            // Heun predictor step
            vcl::kernel_call(predictor_step);

            vcl::queue.finish();

            // update spin fields, external fixed
            vcl::update_spin_fields();

            vcl::queue.finish();

            // Heun corrector step
            vcl::kernel_call(corrector_step);

            vcl::queue.finish();
         }
      }
   }
}
#endif // OPENCL
