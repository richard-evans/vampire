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

namespace vcuda{

#ifdef CUDA

   namespace internal{

      namespace atoms
      {
         /*
          * Thrust array instantiation
          */
         RealArray x_spin_array(0UL);
         RealArray y_spin_array(0UL);
         RealArray z_spin_array(0UL);
         RealArray x_coord_array(0UL);
         RealArray y_coord_array(0UL);
         RealArray z_coord_array(0UL);
         IndexArray type_array(0UL);
         IndexArray cell_array(0UL);
         IndexArray limits(0UL);
         IndexArray neighbours(0UL);
      } /* atoms */

      namespace cells
      {
         RealArray x_coord_array(0UL);
         RealArray y_coord_array(0UL);
         RealArray z_coord_array(0UL);
         RealArray x_mag_array(0UL);
         RealArray y_mag_array(0UL);
         RealArray z_mag_array(0UL);
         RealArray x_field_array(0UL);
         RealArray y_field_array(0UL);
         RealArray z_field_array(0UL);
         RealArray volume_array(0UL);
         IndexArray num_atoms(0UL);
      } /* cells */

      namespace mp
      {
         MaterialParametersArray materials(0UL);
      } /* mp */

      RealArray x_total_spin_field_array(0UL);
      RealArray y_total_spin_field_array(0UL);
      RealArray z_total_spin_field_array(0UL);
      RealArray x_total_external_field_array(0UL);
      RealArray y_total_external_field_array(0UL);
      RealArray z_total_external_field_array(0UL);
      RealArray x_dipolar_field_array(0UL);
      RealArray y_dipolar_field_array(0UL);
      RealArray z_dipolar_field_array(0UL);

      curandState * d_rand_state;

   } // end of internal namespace

#endif

} // end of vcuda namespace
