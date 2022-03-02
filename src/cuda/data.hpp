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

#include "typedefs.hpp"
/*
 * Provide the definition for the material_t class
 */
#include "internal.hpp"
#include "statistics.hpp"


#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda
{
#ifdef CUDA
   namespace internal
   {

      extern cu_real_t *h_x_spin_transfer_buffer;
      extern cu_real_t *h_y_spin_transfer_buffer;
      extern cu_real_t *h_z_spin_transfer_buffer;

      namespace atoms
      {

         extern cu_real_t *d_x_spin;
         extern cu_real_t *d_y_spin;
         extern cu_real_t *d_z_spin;

         extern cu_real_t *d_spin;

         extern cu_real_t *d_x_coord;
         extern cu_real_t *d_y_coord;
         extern cu_real_t *d_z_coord;

         extern int *d_materials;
         extern int *d_cells;

         extern cu_real_t *d_spin_norm;

         extern int *d_limits;
         extern int *d_neighbours;

         /*
         extern cu_real_array_t x_coord_array;
         extern cu_real_array_t y_coord_array;
         extern cu_real_array_t z_coord_array;

         extern cu_index_array_t type_array;

         extern cu_index_array_t cell_array;

         extern cu_index_array_t limits;
         extern cu_index_array_t neighbours;
         */
         /*
          * Unrolled spin norm array
          */

      } /* atoms */


      namespace cells
      {

         extern cu_real_t *d_x_coord;
         extern cu_real_t *d_y_coord;
         extern cu_real_t *d_z_coord;

         extern cu_real_t *d_x_mag;
         extern cu_real_t *d_y_mag;
         extern cu_real_t *d_z_mag;

         extern cu_real_t *d_x_cell_field;
         extern cu_real_t *d_y_cell_field;
         extern cu_real_t *d_z_cell_field;

         extern cu_real_t *d_x_cell_mu0H_field;
         extern cu_real_t *d_y_cell_mu0H_field;
         extern cu_real_t *d_z_cell_mu0H_field;

         extern cu_real_t *d_tensor_xx;
         extern cu_real_t *d_tensor_xy;
         extern cu_real_t *d_tensor_xz;
         extern cu_real_t *d_tensor_yy;
         extern cu_real_t *d_tensor_yz;
         extern cu_real_t *d_tensor_zz;

         extern cu_real_t *d_volume;
         extern cu_real_t *d_num_atoms;

         extern int *d_cell_id_array;
         extern int *d_num_atoms_in_cell;

      } /* cells */

      namespace hamr
      {
         extern cu_real_t d_head_position_x;
         extern cu_real_t d_head_position_y;
         extern cu_real_t d_H_bounds_x;
         extern cu_real_t d_H_bounds_y;
         extern cu_real_t d_laser_sigma_x;
         extern cu_real_t d_laser_sigma_y;
         extern cu_real_t d_NPS;
      } /* hamr */

      namespace mp
      {
         extern material_parameters_t *d_material_params;
      } /* mp */

      extern cu_real_t *d_x_spin_field;
      extern cu_real_t *d_y_spin_field;
      extern cu_real_t *d_z_spin_field;

      extern cu_real_t *d_spin_field;

      extern cu_real_t *d_x_external_field;
      extern cu_real_t *d_y_external_field;
      extern cu_real_t *d_z_external_field;

      extern cu_real_t *d_x_dip_field;
      extern cu_real_t *d_y_dip_field;
      extern cu_real_t *d_z_dip_field;

      extern cu_real_t *d_x_mu0H_dip_field;
      extern cu_real_t *d_y_mu0H_dip_field;
      extern cu_real_t *d_z_mu0H_dip_field;


      /*
       * Required by the total external field calculator
       * and the dipolar field updater
       */

      /*
       * cuRAND states
       */

      extern curandState * d_rand_state;

   } /* internal */
#endif
} /* vcuda */

#endif
