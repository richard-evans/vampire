#include <sstream>
#include <iostream>
#include "atoms.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      namespace rng
      {
         cl::Kernel grng;

         void update_grands(void)
         {
            vcl::kernel_call(grng, vcl::queue, vcl::global, vcl::local);
         }
      }
   }
}

#endif // OPENCL
