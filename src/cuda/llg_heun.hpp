//------------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------
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
         cu_real_t prefactor;
         /**
          * @var lambda * prefactor
          */
         cu_real_t lambda_times_prefactor;
      };

      namespace llg
      {
         /*
          * Private data
          */
         extern bool initialized;

         /*
          * Internal functions
          */
         void __llg_init ();
         void __llg_step ();

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
               cu_real_t dt, size_t num_atoms
               );

         __global__ void llg_heun_corrector_step (
               int * material_id,
               heun_parameters_t * heun_parameters,
               cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
               cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
               cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
               cu_real_t * x_spin_buffer, cu_real_t * y_spin_buffer, cu_real_t * z_spin_buffer,
               cu_real_t * dS_x, cu_real_t * dS_y, cu_real_t * dS_z,
               cu_real_t dt, size_t num_atoms
               );
      } /* llg */

#endif

   } /* internal */
} /* vcuda */

#endif
