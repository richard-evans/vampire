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
      const double total_time = diff.count();
      std::cout << "Simulation took " << total_time << " seconds" << std::endl;

#ifdef OPENCL_TIME_KERNELS
      std::cout << "_Kernel________________|__%_total_time_" << std::endl;
      std::cout << "Spin fields            | " << 100*vcl::time::spin_fields/total_time     << std::endl;
      std::cout << "Matrix multiplication  | " << 100*vcl::time::mat_mul/total_time         << std::endl;
      std::cout << "Random number generator| " << 100*vcl::time::rng/total_time             << std::endl;
      std::cout << "External fields        | " << 100*vcl::time::external_fields/total_time << std::endl;
      std::cout << "Predictor step         | " << 100*vcl::time::predictor_step/total_time  << std::endl;
      std::cout << "Corrector step         | " << 100*vcl::time::corrector_step/total_time  << std::endl;
#endif // OPENCL_TIME_KERNELS

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
