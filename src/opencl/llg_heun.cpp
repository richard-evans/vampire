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
#include "random.hpp"
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

         vcl::Buffer3D spin_buffer_array;

         vcl::Buffer3D dS_array;

         cl::Buffer heun_parameters_device;

         cl::Kernel predictor_step;
         cl::Kernel corrector_step;

         void init(void) noexcept
         {
            const size_t real_buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);
            const size_t num_mats = ::mp::num_materials;

            vcl::llg::spin_buffer_array = vcl::Buffer3D(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);

            vcl::llg::dS_array = vcl::Buffer3D(vcl::context, CL_MEM_READ_WRITE, real_buffer_size);

            vcl::llg::heun_parameters_device = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, num_mats*sizeof(heun_parameter_t));

            std::vector<heun_parameter_t> heun_parameters_host(num_mats);

            for (unsigned i=0; i<num_mats; ++i)
            {
               const double alpha = ::mp::material[i].alpha;
               const double gamma = ::mp::material[i].gamma_rel;
               heun_parameters_host[i].prefactor = -gamma / (1.0 + alpha*alpha);
               heun_parameters_host[i].lambda_times_prefactor = (-gamma * alpha /
                                                                 (1.0 + alpha*alpha));
            }

            vcl::queue.enqueueWriteBuffer(vcl::llg::heun_parameters_device,
                                          CL_FALSE,
                                          0,
                                          num_mats*sizeof(heun_parameter_t),
                                          &heun_parameters_host[0]);

            // Build Kernels
            std::ostringstream opts;
            opts << "-DNUM_ATOMS=" << ::atoms::num_atoms;
            opts << " -DDT=" << ::mp::dt;
            predictor_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                         "llg_heun_predictor_step",
                                                         vcl::context, vcl::default_device, opts.str());

            corrector_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                         "llg_heun_corrector_step",
                                                         vcl::context, vcl::default_device, opts.str());

            vcl::set_kernel_args(predictor_step,
                                 vcl::atoms::type_array,
                                 vcl::llg::heun_parameters_device,
                                 vcl::atoms::spin_array.x(),
                                 vcl::atoms::spin_array.y(),
                                 vcl::atoms::spin_array.z(),
                                 vcl::total_spin_field_array.x(),
                                 vcl::total_spin_field_array.y(),
                                 vcl::total_spin_field_array.z(),
                                 vcl::total_external_field_array.x(),
                                 vcl::total_external_field_array.y(),
                                 vcl::total_external_field_array.z(),
                                 vcl::llg::dS_array.x(),
                                 vcl::llg::dS_array.y(),
                                 vcl::llg::dS_array.z());

            vcl::set_kernel_args(corrector_step,
                                 vcl::atoms::type_array,
                                 vcl::llg::heun_parameters_device,
                                 vcl::atoms::spin_array.x(),
                                 vcl::atoms::spin_array.y(),
                                 vcl::atoms::spin_array.z(),
                                 vcl::total_spin_field_array.x(),
                                 vcl::total_spin_field_array.y(),
                                 vcl::total_spin_field_array.z(),
                                 vcl::total_external_field_array.x(),
                                 vcl::total_external_field_array.y(),
                                 vcl::total_external_field_array.z(),
                                 vcl::llg::spin_buffer_array.x(),
                                 vcl::llg::spin_buffer_array.y(),
                                 vcl::llg::spin_buffer_array.z(),
                                 vcl::llg::dS_array.x(),
                                 vcl::llg::dS_array.y(),
                                 vcl::llg::dS_array.z());

            vcl::llg::initialized = true;
         }

         void step(void) noexcept
         {
            const cl::NDRange global(::atoms::num_atoms);

            const size_t real_buffer_size = ::atoms::num_atoms * sizeof(vcl_real_t);

            vcl::atoms::spin_array.copy_to_dev(vcl::queue, vcl::llg::spin_buffer_array, real_buffer_size);

            // update fields
            vcl::update_spin_fields();
            vcl::rng::update_grands();

            vcl::queue.finish();

            vcl::update_external_fields();

            vcl::queue.finish();

            // Heun predictor step
            vcl::kernel_call(predictor_step, vcl::queue, global, vcl::local);


            vcl::queue.finish();

            // update spin fields, external fixed
            vcl::update_spin_fields();

            vcl::queue.finish();

            // Heun corrector step
            vcl::kernel_call(corrector_step, vcl::queue, global, vcl::local);


            vcl::queue.finish();
         }
      }
   }
}
#endif // OPENCL
