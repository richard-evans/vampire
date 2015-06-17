#ifndef CUDA_INTERNAL_H_
#define CUDA_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// cuda implementation. These functions should
// not be accessed outside of the local temperature pulse code.
//---------------------------------------------------------------------

#include <thrust/copy.h>
#include <thrust/device_vector.h>

/*
 * requesting data strcutures from the main program
 */

#include "../../hdr/atoms.hpp"
#include "../../hdr/cells.hpp"
#include "../../hdr/material.hpp"

namespace cuda{

#ifdef CUDA

   namespace internal{
      struct material_parameters_t {
         double alpha;
         double gamma_rel;
         double mu_s_SI;
         double Klatt_SI;
         double sh2;
         double sh4;
         double sh6;
         double anisotropy_unit_x;
         double anisotropy_unit_y;
         double anisotropy_unit_z;
         double Kc1_SI;
         double temperature;
         double temperature_rescaling_alpha;
         double temperature_rescaling_Tc;
      };

      struct heun_parameters_t {
         /**
          * @var gamma_rel / (1 + alpha ** 2)
          */
         double prefactor;
         /**
          * @var alpha * prefactor
          */
         double lambda_times_prefactor;
      };

      /*
       * Initlialization functions
       */
      bool __initialize_atoms ();
      bool __initialize_fields ();
      bool __initialize_cells ();
      bool __initialize_materials ();
      bool __initialize_topology ();

      /*
       * Shared functors for thrust
       */

      struct plusone_functor
      {
         __host__ __device__
            size_t operator() (const float& item) const
            {
               return item + 1UL;
            }
      };

      /*
       * Shared kernel definitions
       */

      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material, size_t * cell,
            double * x_sp_field, double * y_sp_field, double * z_sp_field
            );

      __global__ void llg_heun_first_kernel (
            double * x_spin, double * y_spin, double * z_spin,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            double dt
            );
     __global__ void llg_heun_scheme(
	    double * x_spin, double * y_spin, double * z_spin,
	    double * x_sp_field, double * y_sp_field, double * z_sp_field,
	    double * x_ext_field, double * y_ext_field, double * z_ext_field,
	    double * x_new_spin, double * y_new_spin, double z_new_spin
            );
#endif
   } // end of iternal namespace
} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
