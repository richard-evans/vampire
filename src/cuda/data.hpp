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
      typedef thrust::device_vector<cu::material_parameters_t>
         MaterialParametersArray;

      namespace atoms
      {
         extern RealArray x_spin_array;
         extern RealArray y_spin_array;
         extern RealArray z_spin_array;

         extern RealArray x_coord_array;
         extern RealArray y_coord_array;
         extern RealArray z_coord_array;

         extern IndexArray type_array;

         extern IndexArray cell_array;

         extern IndexArray limits;
         extern IndexArray neighbours;

      } /* atoms */


      namespace exchange
      {
         extern RealArray Jxx_vals_d;
         extern RealArray Jyy_vals_d;
         extern RealArray Jzz_vals_d;
         /*
          * TODO: Tensor exchanges
          */
      }
      namespace cells
      {
         extern RealArray x_coord_array;
         extern RealArray y_coord_array;
         extern RealArray z_coord_array;

         extern RealArray x_mag_array;
         extern RealArray y_mag_array;
         extern RealArray z_mag_array;

         extern RealArray x_field_array;
         extern RealArray y_field_array;
         extern RealArray z_field_array;

         extern RealArray volume_array;

         extern IndexArray num_atoms;
      } /* cells */

      namespace mp
      {
         extern MaterialParametersArray materials;
      } /* mp */


      extern RealArray x_total_spin_field_array;
      extern RealArray y_total_spin_field_array;
      extern RealArray z_total_spin_field_array;

      extern RealArray x_total_external_field_array;
      extern RealArray y_total_external_field_array;
      extern RealArray z_total_external_field_array;

      /*
       * Required by the total external field calculator
       * and the dipolar field updater
       */
      extern RealArray x_dipolar_field_array;
      extern RealArray y_dipolar_field_array;
      extern RealArray z_dipolar_field_array;

      /*
       * cuRAND states
       */

      extern curandState * d_rand_state;

      namespace llg
      {
         extern RealArray x_spin_prima_array;
         extern RealArray y_spin_prima_array;
         extern RealArray z_spin_prima_array;
      } /* llg */

   } /* internal */
#endif
} /* vcuda */

#endif
