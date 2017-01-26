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
#include "random.hpp"
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
      namespace llg
      {
         cl::Kernel predictor_step;
         cl::Kernel corrector_step;

         void step(void) noexcept
         {
            vcl::atoms::spin_array.copy_to_dev(vcl::queue, vcl::llg::spin_buffer_array);

            // update fields
            vcl::update_spin_fields();
            vcl::rng::update_grands();

            vcl::queue.finish();

            vcl::update_external_fields();

            vcl::queue.finish();

            // Heun predictor step
            vcl::kernel_call(predictor_step, vcl::queue, vcl::global, vcl::local);


            vcl::queue.finish();

            // update spin fields, external fixed
            vcl::update_spin_fields();

            vcl::queue.finish();

            // Heun corrector step
            vcl::kernel_call(corrector_step, vcl::queue, vcl::global, vcl::local);


            vcl::queue.finish();
         }
      }
   }
}
#endif // OPENCL
