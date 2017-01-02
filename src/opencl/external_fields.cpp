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

            const vcl_real_t Hx = sim::H_vec[0] * sim::H_applied;
            const vcl_real_t Hy = sim::H_vec[1] * sim::H_applied;
            const vcl_real_t Hz = sim::H_vec[2] * sim::H_applied;

            vcl::set_kernel_args(update_ext,
                                 vcl::atoms::type_array,
                                 vcl::mp::materials,
                                 vcl::dipolar_field_array.x(),
                                 vcl::dipolar_field_array.y(),
                                 vcl::dipolar_field_array.z(),
                                 vcl::total_external_field_array.x(),
                                 vcl::total_external_field_array.y(),
                                 vcl::total_external_field_array.z(),
                                 vcl::rng::grands,
                                 sim::temperature,
                                 Hx, Hy, Hz);
         }

         const cl::NDRange global(::atoms::num_atoms);

         vcl::kernel_call(update_ext, vcl::queue, global, vcl::local);

         update_dipolar_fields();
      }
   }
}

#endif // OPENCL
