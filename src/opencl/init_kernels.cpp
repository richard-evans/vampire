#include <sstream>
#include <vector>

#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "sim.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "kernels.hpp"
#include "llg_heun.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

static void init_dipole(void)
{
   // dipole calculations enabled?
   if (::sim::hamiltonian_simulation_flags[4]!=1) return;

   std::ostringstream opts;
   opts << "-DN_CELLS=" << ::cells::num_cells;
   opts << " -DN_ATOMS=" << ::atoms::num_atoms;
   vcl::update_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                 "update_dipole_fields",
                                                 vcl::context, vcl::default_device,
                                                 opts.str());
   vcl::update_atm_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                     "update_atm_dipole_fields",
                                                     vcl::context, vcl::default_device,
                                                     opts.str());
   vcl::update_cell_mag = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                      "update_cell_magnetization",
                                                      vcl::context, vcl::default_device,
                                                      opts.str());

   vcl::set_kernel_args(vcl::update_dip,
                        vcl::cells::mag_array.buffer(),
                        vcl::cells::coord_array.buffer(),
                        vcl::cells::volume_array,
                        vcl::cells::field_array.buffer());

   vcl::set_kernel_args(vcl::update_atm_dip,
                        vcl::cells::field_array.buffer(),
                        vcl::dipolar_field_array.buffer(),
                        vcl::atoms::cell_array);

   vcl::set_kernel_args(vcl::update_cell_mag,
                        vcl::atoms::spin_array.buffer(),
                        vcl::atoms::type_array,
                        vcl::atoms::cell_array,
                        vcl::mp::materials,
                        vcl::cells::mag_array.buffer());
}

static void init_external_fields(void)
{
   std::ostringstream opts;
   opts << "-DNUM_ATOMS=" << ::atoms::num_atoms;
   vcl::update_ext = vcl::build_kernel_from_file("src/opencl/cl/external_fields.cl",
                                                 "update_external_fields",
                                                 vcl::context, vcl::default_device,
                                                 opts.str());

   vcl::set_kernel_args(vcl::update_ext,
                        vcl::atoms::type_array,
                        vcl::mp::materials,
                        vcl::dipolar_field_array.buffer(),
                        vcl::total_external_field_array.buffer(),
                        vcl::rng::grands,
                        vcl::real_t(sim::H_vec[0] * sim::H_applied),
                        vcl::real_t(sim::H_vec[1] * sim::H_applied),
                        vcl::real_t(sim::H_vec[2] * sim::H_applied),
                        vcl::real_t(sim::temperature));
}

static void init_exchange(void)
{
   const size_t vsize = ::atoms::neighbour_list_array.size();

   std::ostringstream opts;
   opts << "-DN=" << ::atoms::num_atoms;
   vcl::exchange::calculate_exchange = vcl::build_kernel_from_file("src/opencl/cl/exchange.cl",
                                                                   "calculate_exchange",
                                                                   vcl::context, vcl::default_device,
                                                                   opts.str());
   switch (::atoms::exchange_type)
   {
   case 0:
      // Isotropic
      // Jxx = Jyy = Jzz
      // Jxy = Jxz = Jyx = 0
   {
      std::vector<vcl::real_t> Jxx_vals_h(vsize);

      for (unsigned i=0; i<vsize; ++i)
      {
         const int iid = ::atoms::neighbour_interaction_type_array[i];
         const vcl::real_t Jij = ::atoms::i_exchange_list[iid].Jij;

         Jxx_vals_h[i] = -Jij;
      }

      vcl::exchange::Jxx_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl::real_t));
      vcl::queue.enqueueWriteBuffer(vcl::exchange::Jxx_vals_d, CL_FALSE, 0, vsize*sizeof(vcl::real_t), &Jxx_vals_h[0]);
   }
   break;
   case 1:
      // Vector
      // Jxx != Jyy != Jzz
      // Jxy = Hxz = Jyx = 0
   {
      std::vector<vcl::real_t> Jxx_vals_h(vsize);
      std::vector<vcl::real_t> Jyy_vals_h(vsize);
      std::vector<vcl::real_t> Jzz_vals_h(vsize);

      for (unsigned i=0; i<vsize; ++i)
      {
         const int iid = ::atoms::neighbour_interaction_type_array[i];

         Jxx_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[0];
         Jyy_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[1];
         Jzz_vals_h[i] = - ::atoms::v_exchange_list[iid].Jij[2];
      }

      vcl::exchange::Jxx_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl::real_t));
      vcl::exchange::Jyy_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl::real_t));
      vcl::exchange::Jzz_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, vsize*sizeof(vcl::real_t));

      vcl::queue.enqueueWriteBuffer(vcl::exchange::Jxx_vals_d, CL_FALSE, 0, vsize*sizeof(vcl::real_t), &Jxx_vals_h[0]);
      vcl::queue.enqueueWriteBuffer(vcl::exchange::Jyy_vals_d, CL_FALSE, 0, vsize*sizeof(vcl::real_t), &Jyy_vals_h[0]);
      vcl::queue.enqueueWriteBuffer(vcl::exchange::Jzz_vals_d, CL_FALSE, 0, vsize*sizeof(vcl::real_t), &Jzz_vals_h[0]);
   }
   break;

   case 2:
      // Tensor
      break;
   }

   vcl::set_kernel_args(vcl::exchange::calculate_exchange,
                        ::atoms::exchange_type,
                        vcl::exchange::Jxx_vals_d,
                        vcl::exchange::Jyy_vals_d,
                        vcl::exchange::Jzz_vals_d,
                        vcl::atoms::limits,
                        vcl::atoms::neighbours,
                        vcl::atoms::spin_array.buffer(),
                        vcl::total_spin_field_array.buffer());
}

static void init_rng(void)
{
   std::ostringstream opts;
   opts << "-DN=" << ::atoms::num_atoms*3;
   vcl::rng::grng = vcl::build_kernel_from_file("src/opencl/cl/random.cl",
                                                "gen_grands",
                                                vcl::context, vcl::default_device,
                                                opts.str());

   vcl::set_kernel_args(vcl::rng::grng, vcl::rng::urands, vcl::rng::grands);
}

static void init_llg(void)
{
   const size_t num_mats = ::mp::num_materials;
   const size_t num_atms = ::atoms::num_atoms;

   vcl::llg::spin_buffer_array =
      vcl::Buffer3D<vcl::real_t>(vcl::context, CL_MEM_READ_WRITE, num_atms);

   vcl::llg::dS_array =
      vcl::Buffer3D<vcl::real_t>(vcl::context, CL_MEM_READ_WRITE, num_atms);

   vcl::llg::heun_parameters_device =
      cl::Buffer(vcl::context, CL_MEM_READ_ONLY, num_mats*sizeof(vcl::heun_parameter_t));

   std::vector<vcl::heun_parameter_t> heun_params_host(num_mats);

   for (unsigned i=0; i<num_mats; ++i)
   {
      const double alpha = ::mp::material[i].alpha;
      const double gamma = ::mp::material[i].gamma_rel;

      heun_params_host[i].prefactor = -gamma / (1.0 + alpha*alpha);
      heun_params_host[i].lambda_times_prefactor = (-gamma * alpha /
                                                    (1.0 + alpha*alpha));
   }

   vcl::queue.enqueueWriteBuffer(vcl::llg::heun_parameters_device,
                                 CL_FALSE, 0,
                                 num_mats*sizeof(vcl::heun_parameter_t),
                                 &heun_params_host[0]);

   std::ostringstream opts;
   opts << "-DNUM_ATOMS=" << num_atms;
   opts << " -DDT=" << ::mp::dt;

   vcl::llg::predictor_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                          "llg_heun_predictor_step",
                                                          vcl::context, vcl::default_device,
                                                          opts.str());
   vcl::llg::corrector_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                          "llg_heun_corrector_step",
                                                          vcl::context, vcl::default_device,
                                                          opts.str());

   vcl::set_kernel_args(vcl::llg::predictor_step,
                        vcl::atoms::type_array,
                        vcl::llg::heun_parameters_device,
                        vcl::atoms::spin_array.buffer(),
                        vcl::total_spin_field_array.buffer(),
                        vcl::total_external_field_array.buffer(),
                        vcl::llg::dS_array.buffer());

   vcl::set_kernel_args(vcl::llg::corrector_step,
                        vcl::atoms::type_array,
                        vcl::llg::heun_parameters_device,
                        vcl::atoms::spin_array.buffer(),
                        vcl::total_spin_field_array.buffer(),
                        vcl::total_external_field_array.buffer(),
                        vcl::llg::spin_buffer_array.buffer(),
                        vcl::llg::dS_array.buffer());
}

namespace vopencl
{
   namespace internal
   {
      bool initialize_kernels(void)
      {
         init_dipole();
         init_external_fields();
         init_exchange();
         init_llg();
         init_rng();

         return true;
      }
   }
}

#endif // OPENCL
