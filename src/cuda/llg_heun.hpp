#ifndef CUDA_LLG_HEUN_HPP_
#define CUDA_LLG_HEUN_HPP_

#include "cuda.hpp"
#include "data.hpp"
#include "internal.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda
{
   namespace internal
   {

#ifdef CUDA

      struct heun_parameters_t {
         /**
          * @var gamma_rel / (1 + alpha ** 2)
          */
         double prefactor;
         /**
          * @var lambda * prefactor
          */
         double lambda_times_prefactor;
      };

      typedef thrust::device_vector<heun_parameters_t> HeunParametersArray;

      namespace llg
      {
         /*
          * Private data
          */

         extern bool initialized;
         extern RealArray x_spin_buffer_array;
         extern RealArray y_spin_buffer_array;
         extern RealArray z_spin_buffer_array;

         extern RealArray dS_x_array;
         extern RealArray dS_y_array;
         extern RealArray dS_z_array;

         extern HeunParametersArray heun_parameters;

         /*
          * Internal functions
          */
         void __llg_init ();
         void __llg_step ();

         /*
          * Internal kernels
          */
         __global__ void llg_heun_step (
               double * x_spin, double * y_spin, double * z_spin,
               int * material_id,
               heun_parameters_t * heun_parameters,
               double * x_sp_field, double * y_sp_field, double * z_sp_field,
               double * x_ext_field, double * y_ext_field, double * z_ext_field,
               double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
               double * dS_x, double * dS_y, double * dS_z,
               double dt, size_t num_atoms
               );

         __global__ void llg_heun_scheme (
               double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
               int * material_id,
               heun_parameters_t * heun_parameters,
               double * x_sp_field, double * y_sp_field, double * z_sp_field,
               double * x_ext_field, double * y_ext_field, double * z_ext_field,
               double * x_spin, double * y_spin, double * z_spin,
               double * dS_x, double * dS_y, double * dS_z,
               double dt, size_t num_atoms
               );
      } /* llg */

#endif

   } /* internal */
} /* vcuda */

#endif
