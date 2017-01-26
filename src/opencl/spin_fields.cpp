#include <sstream>

#include "atoms.hpp"

#include "data.hpp"
#include "exchange_fields.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      cl::Kernel update_nexch_spin_fields;

      void update_spin_fields(void) noexcept
      {
         vcl::total_spin_field_array.zero_buffer();

         vcl::kernel_call(update_nexch_spin_fields, vcl::queue, vcl::global, vcl::local);

         vcl::exchange::calculate_exchange_fields();
      }
   }
}

#endif // OPENCL
