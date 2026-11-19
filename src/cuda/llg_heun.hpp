//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015.
//   All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef CUDA_LLG_HEUN_H_
#define CUDA_LLG_HEUN_H_

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
         cu_real_t prefactor;
         /**
          * @var lambda * prefactor
          */
         cu_real_t lambda_times_prefactor;
         /**
          * @var atoms of this material are not integrated (gpu:skip-non-magnetic-atoms)
          */
         int skip;
      };

      namespace llg
      {
         /*
          * Private data
          */
         extern bool initialized;
         extern heun_parameters_t * d_heun_params;
         extern cu_real_t * d_x_spin_buffer;
         extern cu_real_t * d_y_spin_buffer;
         extern cu_real_t * d_z_spin_buffer;
         extern cu_real_t * d_ds_x;
         extern cu_real_t * d_ds_y;
         extern cu_real_t * d_ds_z;

         /*
          * Internal functions
          */
         void __llg_init ();
         void __llg_step ();
         void __llg_finalise ();

         /*
          * Internal kernels
          */
         __global__ void llg_heun_predictor_step (
               int * material_id,
               heun_parameters_t * heun_parameters,
               cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
               cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
               cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
               cu_real_t * dS_x, cu_real_t * dS_y, cu_real_t * dS_z,
               cu_real_t dt, size_t num_atoms, const int * active
               );

         __global__ void llg_heun_corrector_step (
               int * material_id,
               heun_parameters_t * heun_parameters,
               cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
               cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
               cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
               cu_real_t * x_spin_buffer, cu_real_t * y_spin_buffer, cu_real_t * z_spin_buffer,
               cu_real_t * dS_x, cu_real_t * dS_y, cu_real_t * dS_z,
               cu_real_t dt, size_t num_atoms, const int * active
               );
      } /* llg */

#endif

   } /* internal */
} /* vcuda */

#endif
