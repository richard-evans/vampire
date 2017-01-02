#include <vector>
#include <sstream>

#include "atoms.hpp"

#include "internal.hpp"
#include "data.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      namespace exchange
      {
         bool exchange_initialized = false;
         bool J_isot_initialized = false;
         bool J_vect_initialised = false;
         bool J_tens_initialised = false;

         cl::Buffer Jxx_vals_d;
         cl::Buffer Jyy_vals_d;
         cl::Buffer Jzz_vals_d;

         cl::Kernel matmul;

         void initialize_exchange(void) noexcept
         {
            std::vector<vcl_real_t> Jxx_vals_h;
            std::vector<vcl_real_t> Jyy_vals_h;
            std::vector<vcl_real_t> Jzz_vals_h;

            const size_t vsize = ::atoms::neighbour_list_array.size();

            std::ostringstream opts;
            opts << "-DN=" << vsize;
            matmul = vcl::build_kernel_from_file("src/opencl/cl/csrmatmul.cl", "matmul",
                                                 vcl::context, vcl::default_device, opts.str());

            switch(::atoms::exchange_type)
            {
            case 0:
               // Isotropic
               // Jxx = Jyy = Jzz
               // Jxy = Jxz = Jyx = 0
            {

               Jxx_vals_h.resize(vsize, 0);

               for (unsigned i=0; i<vsize; ++i)
               {
                  const int iid = ::atoms::neighbour_interaction_type_array[i];
                  const vcl_real_t Jij = ::atoms::i_exchange_list[iid].Jij;

                  Jxx_vals_h[i] = - Jij;
               }

               Jxx_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl_real_t));
               vcl::queue.enqueueWriteBuffer(Jxx_vals_d, CL_FALSE, 0, vsize*sizeof(vcl_real_t), &Jxx_vals_h[0]);
               vcl::queue.finish();

               J_isot_initialized = true;

               break;
            }

            case 1:
               // Vector
               // Jxx != Jyy != Jzz
               // Jxy = Jxz = Jyx = 0
            {
               Jxx_vals_h.resize(vsize);
               Jyy_vals_h.resize(vsize);
               Jzz_vals_h.resize(vsize);

               for (unsigned i=0; i<vsize; ++i)
               {
                  const int iid = ::atoms::neighbour_interaction_type_array[i];
                  Jxx_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[0];
                  Jyy_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[1];
                  Jzz_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[2];
               }

               Jxx_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl_real_t));
               Jyy_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl_real_t));
               Jzz_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl_real_t));

               vcl::queue.enqueueWriteBuffer(Jxx_vals_d, CL_FALSE, 0, vsize*sizeof(vcl_real_t), &Jxx_vals_h[0]);
               vcl::queue.enqueueWriteBuffer(Jyy_vals_d, CL_FALSE, 0, vsize*sizeof(vcl_real_t), &Jyy_vals_h[0]);
               vcl::queue.enqueueWriteBuffer(Jzz_vals_d, CL_FALSE, 0, vsize*sizeof(vcl_real_t), &Jzz_vals_h[0]);
               vcl::queue.finish();

               J_vect_initialised = true;

               break;
            }

            case 2:
               // Tensor
               break;
            }

            exchange_initialized = true;
         }

         void calculate_exchange_fields(void) noexcept
         {
            if(!exchange_initialized)
               initialize_exchange();

            // rowptrs = vcl::atoms::limits
            // colidxs = vcl::atoms::neighbours

            // convert Jnn from CSR to DIA

            const cl::NDRange global(::atoms::num_atoms);

            switch(::atoms::exchange_type)
            {
            case 0:
               // Isotropic
               // Jxx = Jyy = Jzz
               // Jxy = Jxz = Jyx = 0

               // vcl::x_total_spin_field_array = matmul(Jxx, vcl::atoms::x_spin_array)
               // vcl::y_total_spin_field_array = matmul(Jxx, vcl::atoms::y_spin_array)
               // vcl::z_total_spin_field_array = matmul(Jxx, vcl::atoms::z_spin_array)
               vcl::set_kernel_args(matmul,
                                    Jxx_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.x(),
                                    vcl::total_spin_field_array.x());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);

               vcl::set_kernel_args(matmul,
                                    Jxx_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.y(),
                                    vcl::total_spin_field_array.y());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);


               vcl::set_kernel_args(matmul,
                                    Jxx_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.z(),
                                    vcl::total_spin_field_array.z());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);


               break;
            case 1:
               // Vector
               // Jxx != Jyy != Jzz
               // Jxy = Jxz = Jyx = 0

               // vcl::x_total_field_array = matmul(Jxx, vcl::atoms::x_spin_array)
               // vcl::y_total_field_array = matmul(Jyy, vcl::atoms::y_spin_array)
               // vcl::z_total_field_array = matmul(Jzz, vcl::atoms::z_spin_array)

               vcl::set_kernel_args(matmul,
                                    Jxx_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.x(),
                                    vcl::total_spin_field_array.x());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);


               vcl::set_kernel_args(matmul,
                                    Jyy_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.y(),
                                    vcl::total_spin_field_array.y());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);


               vcl::set_kernel_args(matmul,
                                    Jzz_vals_d,
                                    vcl::atoms::limits,
                                    vcl::atoms::neighbours,
                                    vcl::atoms::spin_array.z(),
                                    vcl::total_spin_field_array.z());
               vcl::kernel_call(matmul, vcl::queue, global, vcl::local);

               break;
            case 2:
               // Tensor
               break;
            }

            vcl::queue.finish();
         }
      }
   }
}

#endif // OPENCL
