//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef CUDA_TYPEDEF_HPP_
#define CUDA_TYPEDEF_HPP_

// C++ standard library headers

// CUDA and thrust headers
#include <curand_kernel.h>
#include "cusparse.h"

// Vampire headers
#include "cuda.hpp"

// Local cuda headers

// Conditional compilation of all cuda code
#ifdef CUDA

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

   // new type definitions
   #ifdef CUDA_DP
      typedef double cu_real_t;
      const cudaDataType CUSPARSE_REAL = CUDA_R_64F;
   #else
      typedef float cu_real_t;
      const cudaDataType CUSPARSE_REAL = CUDA_R_32F;
   #endif

   //typedef cusp::array1d<cu_real_t, cusp::device_memory> cu_real_array_t;
   //typedef cusp::array1d<int, cusp::device_memory> cu_index_array_t;
   //typedef thrust::device_vector<double> cu_real_array_t;
   //typedef thrust::device_vector<int>    cu_index_array_t;

   // Compile-time selectable matrix structure
   /*
   #if CUDA_MATRIX == CSR
      typedef cusp::csr_matrix<int, cu_real_t, cusp::device_memory> cu_exch_mat_t;
   #elif CUDA_MATRIX == DIA
      typedef cusp::dia_matrix<int, cu_real_t, cusp::device_memory> cu_exch_mat_t;
   #elif CUDA_MATRIX == ELL
      typedef cusp::ell_matrix<int, cu_real_t, cusp::device_memory> cu_exch_mat_t;
   #else
      typedef cusp::csr_matrix<int, cu_real_t, cusp::device_memory> cu_exch_mat_t;
   #endif
   */
   typedef cusparseSpMatDescr_t cu_exch_mat_t;

   // struct for material parameters
   struct material_parameters_t {
      cu_real_t alpha;
      cu_real_t gamma_rel;
      cu_real_t mu_s_si;
      cu_real_t i_mu_s_si;
      cu_real_t k_latt;
      cu_real_t sh2;
      cu_real_t sh4;
      cu_real_t sh6;
      cu_real_t anisotropy_unit_x;
      cu_real_t anisotropy_unit_y;
      cu_real_t anisotropy_unit_z;
      cu_real_t applied_field_strength;
      cu_real_t applied_field_unit_x;
      cu_real_t applied_field_unit_y;
      cu_real_t applied_field_unit_z;
      cu_real_t kc4;
      cu_real_t temperature;
      cu_real_t temperature_rescaling_alpha;
      cu_real_t temperature_rescaling_Tc;
      cu_real_t H_th_sigma;
   };

   // Type definition for array of material parameters
   //typedef thrust::device_vector<material_parameters_t> cu_material_array_t;

} // end of internal namespace

} // end of vcuda namespace

#endif // Conditional compilation

#endif // Header guard
