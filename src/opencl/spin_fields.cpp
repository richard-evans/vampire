#include <sstream>

#include "atoms.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      bool compiled_update_spin_fields = false;
      cl::Kernel update_nexch_spin_fields;

      void update_spin_fields()
      {
         cl::CommandQueue write_q(vcl::context, vcl::default_device);
         size_t buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);

         // Zero the field buffers
         vcl_real_t zero = 0.0;
         write_q.enqueueFillBuffer(vcl::total_spin_field_array.x(),
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);
         write_q.enqueueFillBuffer(vcl::total_spin_field_array.y(),
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);
         write_q.enqueueFillBuffer(vcl::total_spin_field_array.z(),
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);

         if (!compiled_update_spin_fields)
         {
            std::ostringstream opts;
            opts << "-DN_ATOMS=" << ::atoms::num_atoms;
            update_nexch_spin_fields = vcl::build_kernel_from_file("cl/spin_fields.cl",
                                                                   "update_nexch_spin_fields",
                                                                   vcl::context, vcl::default_device,
                                                                   opts.str());
            compiled_update_spin_fields = true;
         }

         write_q.finish();

         cl::NDRange global(::atoms::num_atoms);
         cl::NDRange local(0);

         vcl::kernel_call(update_nexch_spin_fields, write_q, global, local,
                          vcl::atoms::type_array,
                          vcl::mp::materials,
                          vcl::atoms::spin_array.x(),
                          vcl::atoms::spin_array.y(),
                          vcl::atoms::spin_array.z(),
                          vcl::total_spin_field_array.x(),
                          vcl::total_spin_field_array.y(),
                          vcl::total_spin_field_array.z());
      }
   }
}

#endif // OPENCL
