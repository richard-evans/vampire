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

#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/tuple.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <curand_kernel.h>

/*
 * requesting data strcutures from the main program
 */

#include "../../hdr/atoms.hpp"
#include "../../hdr/cells.hpp"
#include "../../hdr/demag.hpp"
#include "../../hdr/material.hpp"
#include "../../hdr/sim.hpp"

// Include type definitions for cuda code
#include "typedefs.hpp"

namespace vcuda{

#ifdef CUDA

   namespace internal{

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

      /*
       * Shared functors for thrust
       */

      template <typename T>
         struct plusone_functor
         {
            __host__ __device__
               int operator() (const T& item) const
               {
                  return item + T(1);
               }
         };

      template <typename T>
         struct scalar_product_functor
         : public thrust::binary_function<
           thrust::tuple<T, T, T>, T, thrust::tuple<T, T, T> >
      {

         typedef thrust::tuple<T, T, T> VecT3;

         __host__ __device__
            VecT3 operator () (const VecT3& a, const T& c)
            {
               return thrust::make_tuple(
                     thrust::get<0>(a) * c,
                     thrust::get<1>(a) * c,
                     thrust::get<2>(a) * c);
            }
      };

      template <typename T>
         struct tuple3_plus_functor
         : public thrust::binary_function<
           thrust::tuple<T, T, T>,
           thrust::tuple<T, T, T>,
           thrust::tuple<T, T, T> >
      {
         typedef thrust::tuple<T, T, T> Tuple3;
         __host__ __device__ Tuple3 operator () (
               const Tuple3 & a, const Tuple3 & b)
         {
            return thrust::make_tuple(
                  thrust::get<0>(a) + thrust::get<0>(b),
                  thrust::get<1>(a) + thrust::get<1>(b),
                  thrust::get<2>(a) + thrust::get<2>(b));
         }
      };



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
            int * material,
            material_parameters_t * material_params,
            cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
            cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
            curandState * rand_state,
            cu_real_t global_temperature,
            cu_real_t Hx, cu_real_t Hy, cu_real_t Hz,
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
            cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
            int n_cells
            );

      __global__ void update_atomistic_dipolar_fields (
            cu_real_t * x_cell_field, cu_real_t * y_cell_field, cu_real_t * z_cell_field,
            cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
            int * cells, int n_atoms
            );

      namespace stats{
         extern bool use_cpu;
      }

   } // end of iternal namespace

#endif

} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
