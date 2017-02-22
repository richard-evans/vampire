//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <chrono>
#include <iostream>

#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"

#ifdef OPENCL
namespace vcl = ::vopencl::internal;
#endif

namespace vopencl
{
   void finalize(void)
   {
#ifdef OPENCL

      auto sim_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = sim_end - vcl::time::sim_start;
      std::cout << "Simulation took " << diff.count() << " seconds" << std::endl;

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
