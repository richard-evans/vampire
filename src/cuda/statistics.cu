//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <thrust/sort.h>

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "internal.hpp"
#include "statistics.hpp"

#ifdef CUDA
namespace cu = vcuda::internal;
#endif

namespace vcuda{

   //-------------------------------------------------------------------------------
   // Function to update statistics
   //-------------------------------------------------------------------------------
   void stats_update(){

      #ifdef CUDA


      #endif

      return;
   }

#ifdef CUDA

   namespace internal
   {
      namespace stats
      {

         void __update_stat (
               const IndexArray & mask,
               const RealArray & stat_saturation,
               RealArray& stat,
               RealArray& mean_stat
               )
         {

         }


         __global__ void hist_by_key_small_mask (
               const double * __restrict__ x_spin,
               const double * __restrict__ y_spin,
               const double * __restrict__ z_spin,
               const double * __restrict__ norm_spin,
               const int * __restrict__ mask,
               double * hist,
               int n_bins,
               int n_atoms
               )
         {
            extern __shared__ double block_hist[];

            for (int i = threadIdx.x; i < 4 * n_bins; i += blockDim.x)
            {
               /*
                * Initialize block memory
                */
               block_hist[i] = 0.0;
            }

            __syncthreads ();

            for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
                  i < n_atoms;
                  i += blockDim.x * gridDim.x)
            {
               /*
                * Store stuff in the shared memory
                */
               int bin = mask[i];
               double mu_s = norm_spin[i];
               cu::atomicAdd (block_hist + 4 * bin + 0, x_spin[i] * mu_s);
               cu::atomicAdd (block_hist + 4 * bin + 1, y_spin[i] * mu_s);
               cu::atomicAdd (block_hist + 4 * bin + 2, z_spin[i] * mu_s);
               cu::atomicAdd (block_hist + 4 * bin + 3, mu_s);
            }

            __syncthreads ();

            for (int i = threadIdx.x; i < 4 * n_bins; i += blockDim.x)
            {
               /*
                * Store stuff in the main memory
                */
               cu::atomicAdd (hist + 4 * i + 0, block_hist[4 * i + 0]);
               cu::atomicAdd (hist + 4 * i + 1, block_hist[4 * i + 1]);
               cu::atomicAdd (hist + 4 * i + 2, block_hist[4 * i + 2]);
               cu::atomicAdd (hist + 4 * i + 3, block_hist[4 * i + 3]);
            }

         }


         __global__ void hist_by_key_big_mask (
               const double * __restrict__ x_spin,
               const double * __restrict__ y_spin,
               const double * __restrict__ z_spin,
               const double * __restrict__ norm_spin,
               const int * __restrict__ mask,
               double * hist,
               int n_bins,
               int n_atoms
               )
         {
            for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
                  i < n_atoms;
                  i += blockDim.x * gridDim.x)
            {
               /*
                * Store stuff in the main memory
                */
               int bin = mask[i];
               double mu_s = norm_spin[i];
               cu::atomicAdd (hist + 4 * bin + 0, x_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 1, y_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 2, z_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 3, mu_s);
            }
         }


         __global__ void update_norm_and_accum (
               double * hist,
               double * accum,
               int n_bins
               )
         {
            for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
                  i < n_bins;
                  i += blockDim.x * gridDim.x)
            {
               double mx = hist[4 * i + 0];
               double my = hist[4 * i + 1];
               double mz = hist[4 * i + 2];
               double ms = hist[4 * i + 3];

               double mm = sqrtf (
                     mx * mx +
                     my * my +
                     mz * mz
                     );

               hist[4 * i + 0] = mx / mm;
               hist[4 * i + 1] = my / mm;
               hist[4 * i + 2] = mz / mm;
               hist[4 * i + 3] = mm / ms;

               accum[4 * i + 0] += mx / mm;
               accum[4 * i + 1] += my / mm;
               accum[4 * i + 2] += mz / mm;
               accum[4 * i + 3] += mm / ms;
            }
         }

      } /* stats */
   } /* internal */

#endif

} // end of namespace cuda
