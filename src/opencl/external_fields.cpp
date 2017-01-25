#include <sstream>

#include "atoms.hpp"
#include "sim.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "random.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      cl::Kernel update_ext;

      void update_external_fields(void) noexcept
      {
         const cl::NDRange global(::atoms::num_atoms);

         vcl::kernel_call(update_ext, vcl::queue, global, vcl::local);

         update_dipolar_fields();
      }
   }
}

#endif // OPENCL
