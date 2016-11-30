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
         extern cu_material_array_t materials;
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
