#ifndef CUDA_DATA_HPP_
#define CUDA_DATA_HPP_

#include <thrust/copy.h>
#include <thrust/device_vector.h>

/*
 * Provide the definition for the material_t class
 */
#include "internal.hpp"

namespace vcuda
{
   namespace internal
   {
      typedef thrust::device_vector<double> RealArray;
      typedef thrust::device_vector<size_t> IndexArray;
      typedef thrust::device_vector<material_parameters_t> MaterialParametersArray;

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

      namespace cells
      {
         extern RealArray cell_x_coord_array;
         extern RealArray cell_y_coord_array;
         extern RealArray cell_z_coord_array;

         extern RealArray cell_x_mag_array;
         extern RealArray cell_y_mag_array;
         extern RealArray cell_z_mag_array;

         extern RealArray cell_volume_array;

         extern IndexArray cell_num_atoms;
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

   } /* internal */
} /* vcuda */

#endif
