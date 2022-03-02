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
// not be accessed outside of the cuda code.
//---------------------------------------------------------------------

// Include type definitions for cuda code
#include "typedefs.hpp"

/*
 * requesting data strcutures from the main program
 */

#include "atoms.hpp"
#include "cells.hpp"
#include "demag.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"

#include "cuda_timer.h"

namespace vcuda{

#ifdef CUDA

   namespace internal{

      extern timer cuda_timer;

      /*
       * Thread launch parameters
       */

      extern int block_size;
      extern int grid_size;

      /*
       * Internal data structures
       */

      /*
       * Initlialization functions
       */
      bool __initialize_atoms ();
      bool __initialize_fields ();
      bool __initialize_cells ();
      bool __initialize_dipole ();
      bool __initialize_hamr ();
      bool __initialize_materials ();
      bool __initialize_topology ();
      bool __initialize_curand ();
      bool __initialize_stats ();

      /*
       * Clean up function
       */
      void __finalize ();

      /*
       * Field updates
       */

      void update_spin_fields ();
      void update_external_fields ();
      void update_dipolar_fields ();
      void update_cell_magnetizations ();
      void update_hamr_field ();
      void update_global_thermal_field ();
      void update_applied_fields ();




      /*
       * Shared device functions
       */

      inline __device__ double atomicAdd (double * address, double value)
      {
         unsigned long long int * address_as_ull =
            (unsigned long long int *) address;
         unsigned long long int old = *address_as_ull;
         unsigned long long int assumed;
         do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed,
                  __double_as_longlong(value + __longlong_as_double(assumed)));
         } while (assumed != old);
         return __longlong_as_double(old);
      }

      inline __device__ float atomicAdd (float * address, float value)
      {
          return ::atomicAdd(address, value);
      }

      /*
       * Shared kernel definitions
       */

      __global__ void init_rng (curandState * state, int seed);

      __global__ void update_non_exchange_spin_fields_kernel (
            int * material, material_parameters_t * material_params,
            cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
            cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
            int num_atoms
            );

      __global__ void update_external_fields_kernel (
            cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
            cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
            int num_atoms
            );

      __global__ void update_cell_magnetization (
            cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
            int * material, int * cell,
            material_parameters_t * material_params,
            cu_real_t*  x_mag, cu_real_t * y_mag, cu_real_t * z_mag,
            int num_atoms
            );

      __global__ void update_dipolar_fields (
            cu_real_t * x_mag, cu_real_t * y_mag, cu_real_t * z_mag,
            cu_real_t * x_coord, cu_real_t * y_coord, cu_real_t * z_coord,
            cu_real_t * volume,
            cu_real_t * x_cell_field, cu_real_t * y_cell_field, cu_real_t * z_cell_field,
            cu_real_t * x_cell_mu0H_field, cu_real_t * y_cell_mu0H_field, cu_real_t * z_cell_mu0H_field,
            cu_real_t * d_tensor_xx, cu_real_t * d_tensor_xy, cu_real_t * d_tensor_xz,
            cu_real_t * d_tensor_yy, cu_real_t * d_tensor_yz, cu_real_t * d_tensor_zz,
            int * d_cell_id_array,
            int * d_num_atoms_in_cell,
            int n_local_cells, int n_cells
            );

      __global__ void update_atomistic_dipolar_fields (
            cu_real_t * x_cell_field, cu_real_t * y_cell_field, cu_real_t * z_cell_field,
            cu_real_t * x_cell_mu0H_field, cu_real_t * y_cell_mu0H_field, cu_real_t * z_cell_mu0H_field,
            cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
            cu_real_t * x_mu0H_dip_field, cu_real_t * y_mu0H_dip_field, cu_real_t * z_mu0H_dip_field,
            int * cells, int n_atoms
            );

      __global__ void apply_global_temperature_kernel(
      		cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
      		cu_real_t temperature,
      		curandState * rand_states,
      		material_parameters_t * material_params,
      		int * material,
      		int n_atoms);

      __global__ void update_applied_fields_kernel(
            cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
            const cu_real_t Hx, const cu_real_t Hy, const cu_real_t Hz,
            int *  material, vcuda::internal::material_parameters_t * material_params,
            const int n_atoms);
      

      namespace stats{
         extern bool use_cpu;
      }

   } // end of iternal namespace

#endif

} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
