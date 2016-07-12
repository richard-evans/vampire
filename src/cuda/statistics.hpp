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

         extern cu_index_array_t material_mask;
         extern cu_real_array_t  material_magnetization;
         extern cu_real_array_t  material_mean_magnetization;

         extern cu_index_array_t height_mask;
         extern cu_real_array_t  height_magnetization;
         extern cu_real_array_t  height_mean_magnetization;

         extern cu_index_array_t material_height_mask;
         extern cu_real_array_t  material_height_magnetization;
         extern cu_real_array_t  material_height_mean_magnetization;

         /*
          * Functions required for statistics
          */

         void __update_stat (
               const cu_index_array_t& mask,
               cu_real_array_t& stat,
               cu_real_array_t& mean_stat
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

         __global__ void hist_by_key_small_mask (
               const double * __restrict__ x_spin,
               const double * __restrict__ y_spin,
               const double * __restrict__ z_spin,
               const double * __restrict__ norm_spin,
               const int * __restrict__ mask,
               double * hist,
               int n_bins,
               int n_atoms
               );

         __global__ void hist_by_key_big_mask (
               const double * __restrict__ x_spin,
               const double * __restrict__ y_spin,
               const double * __restrict__ z_spin,
               const double * __restrict__ norm_spin,
               const int * __restrict__ mask,
               double * hist,
               int n_bins,
               int n_atoms
               );

         __global__ void update_norm_and_accum (
               double * hist,
               double * accum,
               int n_bins
               );

      } /* stats */
   } /* internal */
   #endif
} /* vcuda */

#endif
