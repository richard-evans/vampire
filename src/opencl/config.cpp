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
         size_t buff_size = ::atoms::num_atoms * sizeof(vcl::vcl_real_t);
         cl::CommandQueue sync_q(vcl::context, vcl::default_device);

         sync_q.enqueueReadBuffer(vcl::atoms::x_spin_array, CL_FALSE, 0, buff_size, &::atoms::x_spin_array[0]);
         sync_q.enqueueReadBuffer(vcl::atoms::y_spin_array, CL_FALSE, 0, buff_size, &::atoms::y_spin_array[0]);
         sync_q.enqueueReadBuffer(vcl::atoms::z_spin_array, CL_FALSE, 0, buff_size, &::atoms::z_spin_array[0]);

         sync_q.finish();
      }
   }
}

#endif
