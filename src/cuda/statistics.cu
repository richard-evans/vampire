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
               const RealArray& mask,
               const RealArray& stat_saturation,
               RealArray& stat,
               RealArray& mean_stat
               )
         {

            /*
             * Copy the data to the temporary arrays as it will be
             * sorted in place.
             */

            thrust::copy (
                  mask.begin(),
                  mask.end(),
                  temp_mask.begin()
                  );

            thrust::copy (
                  thrust::make_zip_iterator(thrust::make_tuple(
                        cu::atoms::x_spin_array.begin(),
                        cu::atoms::y_spin_array.begin(),
                        cu::atoms::z_spin_array.begin())),
                  thrust::make_zip_iterator(thrust::make_tuple(
                        cu::atoms::x_spin_array.end(),
                        cu::atoms::y_spin_array.end(),
                        cu::atoms::z_spin_array.end())),
                  thrust::make_zip_iterator(thrust::make_tuple(
                        x_temp_spin.begin(),
                        y_temp_spin.begin(),
                        z_temp_spin.begin()))
                  );

            /*
             * Multiply the spins times the spin norm
             */

            thrust::transform(
                  x_temp_spin.begin(),
                  x_temp_spin.end(),
                  spin_norm.begin(),
                  x_temp_spin.begin(),
                  cu::scalar_product_functor<double>()
                  );

            thrust::sort_by_key(
                  temp_mask.begin(),
                  temp_mask.end(),
                  thrust::make_zip_iterator(thrust::make_tuple(
                        x_temp_spin.begin(),
                        y_temp_spin.begin(),
                        z_temp_spin.begin()))
                  );

            size_t mask_size = stat.size() / 4UL;

            thrust::reduce_by_key(
                  temp_mask.begin(),
                  temp_mask.end(),
                  thrust::make_zip_iterator(thrust::make_tuple(
                        x_temp_spin.begin(),
                        y_temp_spin.begin(),
                        z_temp_spin.begin())),
                  thrust::make_discard_iterator(), // We don't need the keys
                  thrust::make_zip_iterator(thrust::make_tuple(
                        stat.begin() + 0UL * mask_size,
                        stat.begin() + 1UL * mask_size,
                        stat.begin() + 2UL * mask_size)),
                  cu::tuple3_plus_functor<double>()
                  );
         }

      } /* stats */
   } /* internal */

#endif

} // end of namespace cuda
