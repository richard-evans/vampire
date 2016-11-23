#include "vopencl.hpp"

#include "data.hpp"

#ifdef OPENCL
namespace vcl = ::vopencl::internal;
#endif

namespace vopencl
{
   void finalize(void)
   {
#ifdef OPENCL

      vcl::atoms::spin_array.free();
      vcl::atoms::coord_array.free();

      vcl::cells::coord_array.free();
      vcl::cells::mag_array.free();
      vcl::cells::field_array.free();

      vcl::total_spin_field_array.free();
      vcl::total_external_field_array.free();
      vcl::dipolar_field_array.free();

#endif //OPENCL      
   }
}
