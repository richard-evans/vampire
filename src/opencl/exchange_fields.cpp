#include <vector>
#include <sstream>

#include "atoms.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      namespace exchange
      {
         cl::Buffer Jxx_vals_d;
         cl::Buffer Jyy_vals_d;
         cl::Buffer Jzz_vals_d;

         cl::Kernel calculate_exchange;

         void calculate_exchange_fields(void) noexcept
         {
            const cl::NDRange global(::atoms::num_atoms);

            vcl::kernel_call(calculate_exchange, vcl::queue, global, vcl::local);

            vcl::queue.finish();
         }
      }
   }
}

#endif // OPENCL
