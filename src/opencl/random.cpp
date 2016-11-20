#include <sstream>

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
         bool compiled_cmwc = false;
         bool compiled_bm = false;

         cl::Kernel cmwc;
         cl::Kernel bm;

         void update_urands(void)
         {
            if (!compiled_cmwc)
            {
               std::ostringstream opts;
               opts << "-DN" << ::atoms::num_atoms*3;
               cmwc = vcl::build_kernel_from_file("cmwc.cl", "cmwc",
                                                  vcl::context, vcl::default_device,
                                                  opts.str());
               compiled_cmwc = true;
            }

            cl::CommandQueue rands_q(vcl::context, vcl::default_device);

            cl::NDRange global(::atoms::num_atoms*3);
            cl::NDRange local(0);

            vcl::kernel_call(cmwc, rands_q, global, local, vcl::rng::urands);
         }

         void update_grands(void)
         {
            if (!compiled_bm)
            {
               std::ostringstream opts;
               opts << "-DN" << ::atoms::num_atoms*3;
               bm = vcl::build_kernel_from_file("boxmuller.cl", "BoxMullerTransform",
                                                vcl::context, vcl::default_device,
                                                opts.str());
               compiled_bm = true;
            }

            update_urands();

            cl::CommandQueue rands_q(vcl::context, vcl::default_device);

            cl::NDRange global(::atoms::num_atoms*3);
            cl::NDRange local(0);

            vcl::kernel_call(bm, rands_q, global, local,
                             vcl::rng::urands,
                             vcl::rng::grands);
         }
      }
   }
}

#endif // OPENCL
