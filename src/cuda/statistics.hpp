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

         extern IndexArray system_mask;
         extern RealArray  system_magnetization;
         extern RealArray  system_mean_magnetization;

         extern IndexArray material_mask;
         extern RealArray  material_magnetization;
         extern RealArray  material_mean_magnetization;

         extern IndexArray height_mask;
         extern RealArray  height_magnetization;
         extern RealArray  height_mean_magnetization;

         extern IndexArray material_height_mask;
         extern RealArray  material_height_magnemagnetization;
         extern RealArray  material_height_mean_mean_magnetization;

         /*
          * Functions required for statistics
          */

         void __update_stat (
               const RealArray& mask,
               RealArray& stat,
               RealArray& mean_stat
               );

         void __get_stat (
               const RealArray& stat,
               const RealArray& mean_stat,
               ::stats::magnetization_statistic_t& local_stat
               );


         void __reset_stat (
               RealArray& stat,
               RealArray& mean_stat
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
