#include "atoms.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace config
   {
      void synchronise(void)
      {
         vcl::atoms::spin_array.copy_to_host(vcl::queue,
                                             ::atoms::x_spin_array,
                                             ::atoms::y_spin_array,
                                             ::atoms::z_spin_array);
      }
   }
}

#endif
