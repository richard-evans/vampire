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
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/tuple.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>

// Vampire headers
#include "cuda.hpp"

// Local cuda headers

// Conditional compilation of all cuda code
#ifdef CUDA

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

   // new type definitions (need to be selectable at compile time)
   typedef double cu_real_t;
   typedef thrust::device_vector<cu_real_t> cu_real_array_t;
   typedef thrust::device_vector<int> cu_index_array_t;

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
      cu_real_t ku;
      cu_real_t anisotropy_unit_x;
      cu_real_t anisotropy_unit_y;
      cu_real_t anisotropy_unit_z;
      cu_real_t applied_field_strength;
      cu_real_t applied_field_unit_x;
      cu_real_t applied_field_unit_y;
      cu_real_t applied_field_unit_z;
      cu_real_t Kc1_SI;
      cu_real_t temperature;
      cu_real_t temperature_rescaling_alpha;
      cu_real_t temperature_rescaling_Tc;
      cu_real_t H_th_sigma;
   };

   // Type definition for array of material parameters
   typedef thrust::device_vector<material_parameters_t> cu_material_array_t;

} // end of internal namespace

} // end of vcuda namespace

#endif // Conditional compilation

#endif // Header guard
