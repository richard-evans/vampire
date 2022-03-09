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

#include "data.hpp"

namespace vcuda{

#ifdef CUDA

   namespace internal{

      cu_real_t *h_x_spin_transfer_buffer;
      cu_real_t *h_y_spin_transfer_buffer;
      cu_real_t *h_z_spin_transfer_buffer;

      namespace atoms
      {
         /*
          * Thrust array instantiation
          */
         cu_real_t *d_x_spin;
         cu_real_t *d_y_spin;
         cu_real_t *d_z_spin;

         cu_real_t *d_spin;

         cu_real_t *d_x_coord;
         cu_real_t *d_y_coord;
         cu_real_t *d_z_coord;

         //cu_real_array_t x_coord_array(0UL);
         //cu_real_array_t y_coord_array(0UL);
         //cu_real_array_t z_coord_array(0UL);

         int *d_materials;
         int *d_cells;
         //cu_index_array_t type_array(0UL);
         //cu_index_array_t cell_array(0UL);

         int *d_limits;
         int *d_neighbours;

         cu_real_t *d_spin_norm;

      } /* atoms */

      namespace cells
      {

         cu_real_t *d_x_coord;
         cu_real_t *d_y_coord;
         cu_real_t *d_z_coord;

         cu_real_t *d_x_mag;
         cu_real_t *d_y_mag;
         cu_real_t *d_z_mag;

         cu_real_t *d_x_cell_field;
         cu_real_t *d_y_cell_field;
         cu_real_t *d_z_cell_field;

         cu_real_t *d_x_cell_mu0H_field;
         cu_real_t *d_y_cell_mu0H_field;
         cu_real_t *d_z_cell_mu0H_field;

         cu_real_t *d_tensor_xx;
         cu_real_t *d_tensor_xy;
         cu_real_t *d_tensor_xz;
         cu_real_t *d_tensor_yy;
         cu_real_t *d_tensor_yz;
         cu_real_t *d_tensor_zz;

         /*
         cu_real_array_t x_coord_array(0UL);
         cu_real_array_t y_coord_array(0UL);
         cu_real_array_t z_coord_array(0UL);

         cu_real_array_t x_mag_array(0UL);
         cu_real_array_t y_mag_array(0UL);
         cu_real_array_t z_mag_array(0UL);

         cu_real_array_t x_field_array(0UL);
         cu_real_array_t y_field_array(0UL);
         cu_real_array_t z_field_array(0UL);
         */

         cu_real_t *d_volume;
         cu_real_t *d_num_atoms;

         //cu_real_array_t volume_array(0UL);
         //cu_index_array_t num_atoms(0UL);
         
         int *d_cell_id_array;
         int *d_num_atoms_in_cell;
 
      } /* cells */

      namespace hamr
      {
         cu_real_t d_head_position_x;
         cu_real_t d_head_position_y;
         cu_real_t d_H_bounds_x;
         cu_real_t d_H_bounds_y;
         cu_real_t d_laser_sigma_x;
         cu_real_t d_laser_sigma_y;
         cu_real_t d_NPS;
      } /* hamr */

      namespace mp
      {
         material_parameters_t *d_material_params;
      } /* mp */

      // Back to the top namespace

      cu_real_t *d_x_spin_field;
      cu_real_t *d_y_spin_field;
      cu_real_t *d_z_spin_field;

      cu_real_t *d_spin_field;

      cu_real_t *d_x_external_field;
      cu_real_t *d_y_external_field;
      cu_real_t *d_z_external_field;

      cu_real_t *d_x_dip_field;
      cu_real_t *d_y_dip_field;
      cu_real_t *d_z_dip_field;

      cu_real_t *d_x_mu0H_dip_field;
      cu_real_t *d_y_mu0H_dip_field;
      cu_real_t *d_z_mu0H_dip_field;


      /*cu_real_array_t x_total_spin_field_array(0UL);
      cu_real_array_t y_total_spin_field_array(0UL);
      cu_real_array_t z_total_spin_field_array(0UL);
      cu_real_array_t x_total_external_field_array(0UL);
      cu_real_array_t y_total_external_field_array(0UL);
      cu_real_array_t z_total_external_field_array(0UL);

      cu_real_array_t x_dipolar_field_array(0UL);
      cu_real_array_t y_dipolar_field_array(0UL);
      cu_real_array_t z_dipolar_field_array(0UL);
      */

      curandState * d_rand_state;

      namespace stats
      {

         bool use_cpu = false;
         long counter(0L);

      } /* stats */

   } // end of internal namespace

#endif

} // end of vcuda namespace
