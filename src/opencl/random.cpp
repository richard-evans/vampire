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
            // launch N/2 work items as each gens 2 numbers
            const cl::NDRange global((::atoms::num_atoms*3)/2);

            vcl::kernel_call(grng, vcl::queue, global, vcl::local);
         }
      }
   }
}

#endif // OPENCL
