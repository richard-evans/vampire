#include <sstream>

#include "atoms.hpp"

#include "data.hpp"
#include "exchange_fields.hpp"
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

      void update_spin_fields(void) noexcept
      {
         vcl::total_spin_field_array.zero_buffer();

         if (!compiled_update_spin_fields)
         {
            std::ostringstream opts;
            opts << "-DNUM_ATOMS=" << ::atoms::num_atoms;
            update_nexch_spin_fields = vcl::build_kernel_from_file("src/opencl/cl/spin_fields.cl",
                                                                   "update_nexch_spin_fields",
                                                                   vcl::context, vcl::default_device,
                                                                   opts.str());
            compiled_update_spin_fields = true;

            vcl::set_kernel_args(update_nexch_spin_fields,
                                 vcl::atoms::type_array,
                                 vcl::mp::materials,
                                 vcl::atoms::spin_array.buffer(),
                                 vcl::total_spin_field_array.buffer());
         }

         cl::NDRange global(::atoms::num_atoms);

         vcl::kernel_call(update_nexch_spin_fields, vcl::queue, global, vcl::local);


         vcl::exchange::calculate_exchange_fields();
      }
   }
}

#endif // OPENCL
