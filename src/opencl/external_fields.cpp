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

      void update_external_fields(void)
      {
         if (!compiled_update_external_fields)
         {
            update_ext = vcl::build_kernel_from_file("cl/external_fields.cl",
                                                     "update_external_fields",
                                                     vcl::context, vcl::default_device);
            compiled_update_external_fields = true;
         }

         cl::CommandQueue update_q(vcl::context, vcl::default_device);

         cl::NDRange global(::atoms::num_atoms);
         cl::NDRange local(0);

         vcl::rng::update_grands();

         vcl_real_t Hx = sim::H_vec[0] * sim::H_applied;
         vcl_real_t Hy = sim::H_vec[1] * sim::H_applied;
         vcl_real_t Hz = sim::H_vec[2] * sim::H_applied;

         vcl::kernel_call(update_ext, update_q, global, local,
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

         update_dipolar_fields();
      }
   }
}

#endif // OPENCL
