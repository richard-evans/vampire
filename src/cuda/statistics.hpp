#ifndef CUDA_STATISTICS_HPP_
#define CUDA_STATISTICS_HPP_

#include "cuda.hpp"
#include "internal.hpp"
#include "data.hpp"

#include "stats.hpp"

namespace vcuda
{
   #ifdef CUDA
   namespace internal
   {
      namespace stats
      {
         /*
          * Arrays required for statistics
          */

         extern long counter;

         extern cu_index_array_t system_mask;
         extern cu_real_array_t  system_magnetization;
         extern cu_real_array_t  system_mean_magnetization;
         extern int system_mask_size;

         extern cu_index_array_t material_mask;
         extern cu_real_array_t  material_magnetization;
         extern cu_real_array_t  material_mean_magnetization;
         extern int material_mask_size;

         extern cu_index_array_t height_mask;
         extern cu_real_array_t  height_magnetization;
         extern cu_real_array_t  height_mean_magnetization;
         extern int height_mask_size;

         extern cu_index_array_t material_height_mask;
         extern cu_real_array_t  material_height_magnetization;
         extern cu_real_array_t  material_height_mean_magnetization;
         extern int material_height_mask_size;

         /*
          * Functions required for statistics
          */

         void __update_stat (
               const cu_index_array_t& mask,
               cu_real_array_t& stat,
               cu_real_array_t& mean_stat,
               int mask_size
               );

         void __get_stat (
               const cu_real_array_t& stat,
               const cu_real_array_t& mean_stat,
               ::stats::magnetization_statistic_t& local_stat
               );


         void __reset_stat (
               cu_real_array_t& stat,
               cu_real_array_t& mean_stat
               );

         /*
          * Kerenels required
          */

         /**
          * This one is super cool but a bit fragile, you should set
          * grid size according to the number of bins in the histogram.
          * Each block will be responsible for updating one of the fields
          * of one of the classes in the histogram.
          */
         __global__ void hist_by_key_smaller_mask (
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               );

         __global__ void hist_by_key_small_mask (
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               );

         __global__ void hist_by_key_big_mask (
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               );

         __global__ void update_norm_and_accum (
               cu_real_t * hist,
               cu_real_t * accum,
               int n_bins
               );

      } /* stats */
   } /* internal */
   #endif
} /* vcuda */

#endif
