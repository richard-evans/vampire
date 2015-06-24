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

#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

/*
 * requesting data strcutures from the main program
 */

#include "../../hdr/atoms.hpp"
#include "../../hdr/cells.hpp"
#include "../../hdr/material.hpp"

namespace vcuda{

#ifdef CUDA

   namespace internal{

      typedef double RealType;

      /*
       * Thread launch parameters
       */

      extern size_t block_size;
      extern size_t grid_size;

      /*
       * Internal data structures
       */

      struct material_parameters_t {
         double alpha;
         double gamma_rel;
         double mu_s_SI;
         double Klatt_SI;
         double sh2;
         double sh4;
         double sh6;
         double ku;
         double ku2;
         double ku3;
         double anisotropy_unit_x;
         double anisotropy_unit_y;
         double anisotropy_unit_z;
         double applied_field_strength;
         double applied_field_unit_x;
         double applied_field_unit_y;
         double applied_field_unit_z;
         double Kc1_SI;
         double temperature;
         double temperature_rescaling_alpha;
         double temperature_rescaling_Tc;
         double H_th_sigma;
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
      bool __initialize_curand ();

      /*
       * Clean up function
       */
      void __finalize ();


      /*
       * Field updates
       */

      void update_spin_fields ();
      void update_external_fields ();

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
       * Shared device functions
       */

      __device__ double atomicAdd (double * address, double value);

      /*
       * Shared kernel definitions
       */

      __global__ void init_rng (curandState * state, size_t seed);

      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material, material_parameters_t * material_params,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            size_t num_atoms
            );

      __global__ void update_external_fields (
            size_t * material, size_t * cell,
            material_parameters_t * material_params,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            curandState * rand_state,
            size_t num_atoms
            );

      __global__ void update_cell_magnetization (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material, size_t * cell,
            material_parameters_t * material_params,
            double * x_mag, double * y_mag, double * z_mag,
            size_t num_atoms
            );

      __global__ void update_dipolar_fields (
            double * x_mag, double * y_mag, double * z_mag,
            double * x_coord, double * y_coord, double * z_coord,
            double * volume, double prefactor, size_t n_cells
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
   } // end of iternal namespace

#endif

} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
