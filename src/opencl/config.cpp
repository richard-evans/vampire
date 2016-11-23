#include "atoms.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace config
   {
      void synchronise(void)
      {
         cl::CommandQueue sync_q(vcl::context, vcl::default_device);

         vcl::atoms::spin_array.copy_to_host(sync_q,
                                             ::atoms::x_spin_array,
                                             ::atoms::y_spin_array,
                                             ::atoms::z_spin_array);

         sync_q.finish();
      }
   }
}

#endif
