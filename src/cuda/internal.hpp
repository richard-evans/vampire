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
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/tuple.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>

/*
 * requesting data strcutures from the main program
 */

#include "../../hdr/atoms.hpp"
#include "../../hdr/cells.hpp"
#include "../../hdr/demag.hpp"
#include "../../hdr/material.hpp"
#include "../../hdr/sim.hpp"

namespace vcuda{

#ifdef CUDA

   namespace internal{

      typedef double RealType;

      /*
       * Thread launch parameters
       */

      extern int block_size;
      extern int grid_size;

      /*
       * Internal data structures
       */

      struct material_parameters_t {
         double alpha;
         double gamma_rel;
         double mu_s_si;
         double i_mu_s_si;
         double k_latt;
         double sh2;
         double sh4;
         double sh6;
         double ku;
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

      __device__ double atomicAdd (double * address, double value);

      /*
       * Shared kernel definitions
       */

      __global__ void init_rng (curandState * state, int seed);

      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            int * material, material_parameters_t * material_params,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            int num_atoms
            );

      __global__ void update_external_fields (
            int * material, int * cell,
            material_parameters_t * material_params,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            curandState * rand_state,
            int num_atoms
            );

      __global__ void update_cell_magnetization (
            double * x_spin, double * y_spin, double * z_spin,
            int * material, int * cell,
            material_parameters_t * material_params,
            double * x_mag, double * y_mag, double * z_mag,
            int num_atoms
            );

      __global__ void update_dipolar_fields (
            double * x_mag, double * y_mag, double * z_mag,
            double * x_coord, double * y_coord, double * z_coord,
            double * volume,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            int n_cells
            );

      __global__ void update_atomistic_dipolar_fields (
            double * x_cell_field, double * y_cell_field, double * z_cell_field,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            int * cells, int n_atoms
            );

   } // end of iternal namespace

#endif

} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
