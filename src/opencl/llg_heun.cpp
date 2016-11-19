#include <vector>
#include <sstream>

#include "atoms.hpp"
#include "material.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "llg_heun.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL
namespace vcl = ::vopencl::internal;
#endif

namespace vopencl
{
   void llg_heun(void)
   {
#ifdef OPENCL
      if (!vcl::llg::initialized)
         vcl::llg::init();

      vcl::llg::step();
#endif // OPENCL
   }
}

#ifdef OPENCL
namespace vopencl
{
   namespace internal
   {
      namespace llg
      {
         bool initialized = false;

         cl::Buffer x_spin_array;
         cl::Buffer y_spin_array;
         cl::Buffer z_spin_array;
         cl::Buffer dS_x_array;
         cl::Buffer dS_y_array;
         cl::Buffer dS_z_array;
         cl::Buffer heun_parameters_device;

         cl::Kernel predictor_step;
         cl::Kernel corrector_step;

         void init(void)
         {
            size_t real_buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);
            size_t num_mats = ::mp::num_materials;

            vcl::llg::x_spin_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);
            vcl::llg::y_spin_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);
            vcl::llg::z_spin_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);

            vcl::llg::dS_x_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);
            vcl::llg::dS_y_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);
            vcl::llg::dS_z_array = cl::Buffer(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);

            vcl::llg::heun_parameters_device = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, num_mats*sizeof(heun_parameter_t));

            std::vector<heun_parameter_t> heun_parameters_host(num_mats);

            for (unsigned i=0; i<num_mats; ++i)
            {
               double alpha = ::mp::material[i].alpha;
               double gamma = ::mp::material[i].gamma_rel;
               heun_parameters_host[i].prefactor = -gamma / (1.0 + alpha*alpha);
               heun_parameters_host[i].lambda_times_prefactor = -gamma * alpha / (1.0 + alpha*alpha);
            }

            cl::CommandQueue write_q(vcl::context, vcl::default_device);
            write_q.enqueueWriteBuffer(vcl::llg::heun_parameters_device, CL_TRUE, 0, num_mats*sizeof(heun_parameter_t), &heun_parameters_host[0]);
            write_q.finish();

            // Build Kernels
            std::ostringstream opts;
            opts << "-DNUM_ATOMS=" << ::atoms::num_atoms;
            predictor_step = vcl::build_kernel_from_file("llg_heun.cl", "llg_heun_predictor_step", vcl::context, vcl::default_device, opts.str());
            corrector_step = vcl::build_kernel_from_file("llg_heun.cl", "llg_heun_corrector_step", vcl::context, vcl::default_device, opts.str());

            vcl::llg::initialized = true;
         }

         void step(void)
         {
            size_t real_buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);
            cl::CommandQueue write_q(vcl::context, vcl::default_device);
            write_q.enqueueCopyBuffer(vcl::atoms::x_spin_array, vcl::llg::x_spin_array, 0, 0, real_buffer_size);
            write_q.enqueueCopyBuffer(vcl::atoms::y_spin_array, vcl::llg::y_spin_array, 0, 0, real_buffer_size);
            write_q.enqueueCopyBuffer(vcl::atoms::z_spin_array, vcl::llg::z_spin_array, 0, 0, real_buffer_size);
         }
      }
   }
}
#endif // OPENCL
