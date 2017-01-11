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
      bool compiled_update_external_fields = false;
      cl::Kernel update_ext;

      void update_external_fields(void) noexcept
      {
         if (!compiled_update_external_fields)
         {
            std::ostringstream opts;
            opts << "-DNUM_ATOMS=" << ::atoms::num_atoms;
            update_ext = vcl::build_kernel_from_file("src/opencl/cl/external_fields.cl",
                                                     "update_external_fields",
                                                     vcl::context, vcl::default_device, opts.str());
            compiled_update_external_fields = true;

            // pack temperature and applied field into one object for minimum global reads
#ifdef OPENCL_DP
            const cl_double4 sys_params =
#else
            const cl_float4 sys_params =
#endif // OPENCL_DP
               {
                  vcl_real_t(sim::H_vec[0] * sim::H_applied),
                  vcl_real_t(sim::H_vec[1] * sim::H_applied),
                  vcl_real_t(sim::H_vec[2] * sim::H_applied),
                  vcl_real_t(sim::temperature)
               };

            vcl::set_kernel_args(update_ext,
                                 vcl::atoms::type_array,
                                 vcl::mp::materials,
                                 vcl::dipolar_field_array.buffer(),
                                 vcl::total_external_field_array.buffer(),
                                 vcl::rng::grands,
                                 sys_params);
         }

         const cl::NDRange global(::atoms::num_atoms);

         vcl::kernel_call(update_ext, vcl::queue, global, vcl::local);

         update_dipolar_fields();
      }
   }
}

#endif // OPENCL
