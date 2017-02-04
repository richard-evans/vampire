//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <vector>

#include "atoms.hpp"
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
   if (::sim::hamiltonian_simulation_flags[4]==1)
   {
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
}

static void init_external_fields(void)
{
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

   switch (::atoms::exchange_type)
   {
   case 0:
      // Isotropic
      // Jxx = Jyy = Jzz
      // Jxy = Jxz = Jyx = 0
   {

#ifdef USE_VECTOR_TYPE
      std::vector<vcl::real_t3> J_vals_h(vsize);
#else
      std::vector<vcl::real_t> J_vals_h(3*vsize);
#endif

      for (unsigned i=0; i<vsize; ++i)
      {
         const int iid = ::atoms::neighbour_interaction_type_array[i];
         const vcl::real_t Jij = ::atoms::i_exchange_list[iid].Jij;

#ifdef USE_VECTOR_TYPE
         J_vals_h[i] = vcl::real_t3{-Jij, -Jij, -Jij};
#else
         const unsigned xxi = 3*i+0;
         const unsigned yyi = 3*i+1;
         const unsigned zzi = 3*i+2;
         J_vals_h[xxi] = J_vals_h[yyi] = J_vals_h[zzi] = - Jij;
#endif
      }

      size_t buffer_size = J_vals_h.size() * sizeof(J_vals_h[0]);
      vcl::exchange::J_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, buffer_size);
      vcl::queue.enqueueWriteBuffer(vcl::exchange::J_vals_d, CL_FALSE, 0, buffer_size, J_vals_h.data());
   }
   break;
   case 1:
      // Vector
      // Jxx != Jyy != Jzz
      // Jxy = Hxz = Jyx = 0
   {
#ifdef USE_VECTOR_TYPE
      std::vector<vcl::real_t3> J_vals_h(vsize);
#else
      std::vector<vcl::real_t> J_vals_h(3*vsize);
#endif

      for (unsigned i=0; i<vsize; ++i)
      {
         const int iid = ::atoms::neighbour_interaction_type_array[i];

         const auto &vel = ::atoms::v_exchange_list[iid];

#ifdef USE_VECTOR_TYPE
         J_vals_h[i] = {-vel.Jij[0], -vel.Jij[1], -vel.Jij[2]};
#else
         const unsigned xxi = 3*i+0;
         const unsigned yyi = 3*i+1;
         const unsigned zzi = 3*i+2;

         J_vals_h[xxi] = - vel.Jij[0];
         J_vals_h[yyi] = - vel.Jij[1];
         J_vals_h[zzi] = - vel.Jij[2];
#endif
      }

      const size_t sz = J_vals_h.size() * sizeof J_vals_h[0];
      vcl::exchange::J_vals_d = cl::Buffer(vcl::context, CL_MEM_READ_ONLY, sz);

      vcl::queue.enqueueWriteBuffer(vcl::exchange::J_vals_d, CL_FALSE, 0, sz, J_vals_h.data());
   }
   break;

   case 2:
      // Tensor
      break;
   }

   vcl::set_kernel_args(vcl::exchange::calculate_exchange,
                        vcl::exchange::J_vals_d,
                        vcl::atoms::limits,
                        vcl::atoms::neighbours,
                        vcl::atoms::spin_array.buffer(),
                        vcl::total_spin_field_array.buffer(),
                        ::atoms::num_atoms);
}

static void init_spin_fields(void)
{
   vcl::set_kernel_args(vcl::update_nexch_spin_fields,
                        vcl::atoms::type_array,
                        vcl::mp::materials,
                        vcl::atoms::spin_array.buffer(),
                        vcl::total_spin_field_array.buffer());
}

static void init_rng(void)
{
   vcl::set_kernel_args(vcl::rng::grng, vcl::rng::state, vcl::rng::grands);
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
         init_spin_fields();
         init_llg();
         init_rng();

         return true;
      }
   }
}

#endif // OPENCL
