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
            size_t buff_size = ::atoms::num_atoms * sizeof(vcl::vcl_real_t);

            cl::CommandQueue stats_q(vcl::context, vcl::default_device);

            stats_q.enqueueReadBuffer(vcl::atoms::x_spin_array, CL_FALSE, 0, buff_size, &::atoms::x_spin_array[0]);
            stats_q.enqueueReadBuffer(vcl::atoms::y_spin_array, CL_FALSE, 0, buff_size, &::atoms::y_spin_array[0]);
            stats_q.enqueueReadBuffer(vcl::atoms::z_spin_array, CL_FALSE, 0, buff_size, &::atoms::z_spin_array[0]);

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
