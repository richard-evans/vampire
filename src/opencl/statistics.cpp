//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "atoms.hpp"
#include "stats.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "statistics.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace stats
   {
      void update(void)
      {
         if (vcl::stats::use_cpu)
         {
            const cl::CommandQueue stats_q(vcl::context, vcl::default_device);

            vcl::atoms::spin_array.copy_to_host(stats_q,
                                                ::atoms::x_spin_array,
                                                ::atoms::y_spin_array,
                                                ::atoms::z_spin_array);

            if (::stats::calculate_system_magnetization)
               ::stats::system_magnetization.calculate_magnetization(::atoms::x_spin_array,
                                                                     ::atoms::y_spin_array,
                                                                     ::atoms::z_spin_array,
                                                                     ::atoms::m_spin_array);
            if (::stats::calculate_material_magnetization)
               ::stats::material_magnetization.calculate_magnetization(::atoms::x_spin_array,
                                                                       ::atoms::y_spin_array,
                                                                       ::atoms::z_spin_array,
                                                                       ::atoms::m_spin_array);
            if (::stats::calculate_height_magnetization)
               ::stats::height_magnetization.calculate_magnetization(::atoms::x_spin_array,
                                                                     ::atoms::y_spin_array,
                                                                     ::atoms::z_spin_array,
                                                                     ::atoms::m_spin_array);
            if (::stats::calculate_material_height_magnetization)
               ::stats::material_height_magnetization.calculate_magnetization(::atoms::x_spin_array,
                                                                              ::atoms::y_spin_array,
                                                                              ::atoms::z_spin_array,
                                                                              ::atoms::m_spin_array);
            return; // Don't do OpenCL version
         }
      }

      void get(void)
      {
         if (vcl::stats::use_cpu) return;
      }

      void reset(void)
      {
         if (vcl::stats::use_cpu)
         {
            if (::stats::calculate_system_magnetization)
               ::stats::system_magnetization.reset_magnetization_averages();

            if (::stats::calculate_material_magnetization)
               ::stats::material_magnetization.reset_magnetization_averages();

            if (::stats::calculate_height_magnetization)
               ::stats::height_magnetization.reset_magnetization_averages();

            if (::stats::calculate_material_height_magnetization)
               ::stats::material_height_magnetization.reset_magnetization_averages();

            return;
         }
      }
   }
}
#endif // OPENCL
