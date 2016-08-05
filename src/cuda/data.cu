//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

/**
 * @brief This file provides object definitions for the global
 *        arrays required by the CUDA implementation.
 */

#include "internal.hpp"
#include "data.hpp"
#include "statistics.hpp"

namespace vcuda{

#ifdef CUDA

   namespace internal{

      namespace atoms
      {
         /*
          * Thrust array instantiation
          */
         cu_real_array_t x_spin_array(0UL);
         cu_real_array_t y_spin_array(0UL);
         cu_real_array_t z_spin_array(0UL);
         cu_real_array_t x_coord_array(0UL);
         cu_real_array_t y_coord_array(0UL);
         cu_real_array_t z_coord_array(0UL);
         cu_index_array_t type_array(0UL);
         cu_index_array_t cell_array(0UL);
         cu_index_array_t limits(0UL);
         cu_index_array_t neighbours(0UL);

         cu_real_array_t spin_norm_array(0UL);

      } /* atoms */

      namespace exchange
      {
         cu_real_array_t Jxx_vals_d(0UL);
         cu_real_array_t Jyy_vals_d(0UL);
         cu_real_array_t Jzz_vals_d(0UL);

         cu_exch_mat_t J_xx_mat_d;
         cu_exch_mat_t J_yy_mat_d;
         cu_exch_mat_t J_zz_mat_d;
      }

      namespace cells
      {
         cu_real_array_t x_coord_array(0UL);
         cu_real_array_t y_coord_array(0UL);
         cu_real_array_t z_coord_array(0UL);
         cu_real_array_t x_mag_array(0UL);
         cu_real_array_t y_mag_array(0UL);
         cu_real_array_t z_mag_array(0UL);
         cu_real_array_t x_field_array(0UL);
         cu_real_array_t y_field_array(0UL);
         cu_real_array_t z_field_array(0UL);
         cu_real_array_t volume_array(0UL);
         cu_index_array_t num_atoms(0UL);
      } /* cells */

      namespace mp
      {
         cu_material_array_t materials(0UL);
      } /* mp */

      cu_real_array_t x_total_spin_field_array(0UL);
      cu_real_array_t y_total_spin_field_array(0UL);
      cu_real_array_t z_total_spin_field_array(0UL);
      cu_real_array_t x_total_external_field_array(0UL);
      cu_real_array_t y_total_external_field_array(0UL);
      cu_real_array_t z_total_external_field_array(0UL);
      cu_real_array_t x_dipolar_field_array(0UL);
      cu_real_array_t y_dipolar_field_array(0UL);
      cu_real_array_t z_dipolar_field_array(0UL);

      curandState * d_rand_state;

      namespace stats
      {

         bool use_cpu = false;
         long counter(0L);

         cu_index_array_t system_mask(0UL);
         cu_real_array_t  system_magnetization(0UL);
         cu_real_array_t  system_mean_magnetization(0UL);
         int system_mask_size(0UL);
         cu_index_array_t material_mask(0UL);
         cu_real_array_t  material_magnetization(0UL);
         cu_real_array_t  material_mean_magnetization(0UL);
         int material_mask_size(0UL);
         cu_index_array_t height_mask(0UL);
         cu_real_array_t  height_magnetization(0UL);
         cu_real_array_t  height_mean_magnetization(0UL);
         int height_mask_size(0UL);
         cu_index_array_t material_height_mask(0UL);
         cu_real_array_t  material_height_magnetization(0UL);
         cu_real_array_t  material_height_mean_magnetization(0UL);
         int material_height_mask_size(0UL);
      } /* stats */

   } // end of internal namespace

#endif

} // end of vcuda namespace
