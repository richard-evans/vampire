//------------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef CUDA_DATA_HPP_
#define CUDA_DATA_HPP_

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

/*
 * Provide the definition for the material_t class
 */
#include "internal.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda
{
#ifdef CUDA
   namespace internal
   {
      typedef thrust::device_vector<double> RealArray;
      typedef thrust::device_vector<int> IndexArray;
      typedef thrust::device_vector<cu::material_parameters_t> MaterialParametersArray;

      // new type definitions (need to be selectable at compile time)
      typedef double cu_real_t;
      typedef thrust::device_vector<double> cu_real_array_t;
      typedef thrust::device_vector<int> cu_index_array_t;
      typedef thrust::device_vector<cu::material_parameters_t> cu_material_array_t;

      namespace atoms
      {
         extern cu_real_array_t x_spin_array;
         extern cu_real_array_t y_spin_array;
         extern cu_real_array_t z_spin_array;

         extern cu_real_array_t x_coord_array;
         extern cu_real_array_t y_coord_array;
         extern cu_real_array_t z_coord_array;

         extern cu_index_array_t type_array;

         extern cu_index_array_t cell_array;

         extern cu_index_array_t limits;
         extern cu_index_array_t neighbours;

         /*
          * Unrolled spin norm array
          */
         extern cu_real_array_t spin_norm_array;

      } /* atoms */


      namespace exchange
      {
         extern cu_real_array_t Jxx_vals_d;
         extern cu_real_array_t Jyy_vals_d;
         extern cu_real_array_t Jzz_vals_d;
         /*
          * TODO: Tensor exchanges
          */
      }
      namespace cells
      {
         extern cu_real_array_t x_coord_array;
         extern cu_real_array_t y_coord_array;
         extern cu_real_array_t z_coord_array;

         extern cu_real_array_t x_mag_array;
         extern cu_real_array_t y_mag_array;
         extern cu_real_array_t z_mag_array;

         extern cu_real_array_t x_field_array;
         extern cu_real_array_t y_field_array;
         extern cu_real_array_t z_field_array;

         extern cu_real_array_t volume_array;

         extern cu_index_array_t num_atoms;
      } /* cells */

      namespace mp
      {
         extern MaterialParametersArray materials;
      } /* mp */


      extern cu_real_array_t x_total_spin_field_array;
      extern cu_real_array_t y_total_spin_field_array;
      extern cu_real_array_t z_total_spin_field_array;

      extern cu_real_array_t x_total_external_field_array;
      extern cu_real_array_t y_total_external_field_array;
      extern cu_real_array_t z_total_external_field_array;

      /*
       * Required by the total external field calculator
       * and the dipolar field updater
       */
      extern cu_real_array_t x_dipolar_field_array;
      extern cu_real_array_t y_dipolar_field_array;
      extern cu_real_array_t z_dipolar_field_array;

      /*
       * cuRAND states
       */

      extern curandState * d_rand_state;

   } /* internal */
#endif
} /* vcuda */

#endif
