#include "atoms.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      void update_spin_fields()
      {
         cl::CommandQueue write_q(vcl::context, vcl::default_device);
         size_t buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);

         // Zero the field buffers
         vcl_real_t zero = 0.0;
         write_q.enqueueFillBuffer(vcl::x_total_spin_field_array,
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);
         write_q.enqueueFillBuffer(vcl::y_total_spin_field_array,
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);
         write_q.enqueueFillBuffer(vcl::z_total_spin_field_array,
                                   &zero,
                                   sizeof(vcl_real_t),
                                   buffer_size);
      }
   }
}

#endif // OPENCL
