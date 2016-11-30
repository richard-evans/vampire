//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "internal.hpp"
#include "statistics.hpp"

#ifdef CUDA
namespace cu = vcuda::internal;
#endif

namespace vcuda{

#ifdef CUDA

   namespace stats
   {
         void update ()
         {

            // If enabled use CPU to calculate statistics by copying data from GPU
            if(vcuda::internal::stats::use_cpu){

				   // copy spin data to CPU
               thrust::copy(internal::atoms::x_spin_array.begin(),internal::atoms::x_spin_array.end(),::atoms::x_spin_array.begin());
               thrust::copy(internal::atoms::y_spin_array.begin(),internal::atoms::y_spin_array.end(),::atoms::y_spin_array.begin());
               thrust::copy(internal::atoms::z_spin_array.begin(),internal::atoms::z_spin_array.end(),::atoms::z_spin_array.begin());

               // call cpu statistics functions
               if(::stats::calculate_system_magnetization)          ::stats::system_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_material_magnetization)        ::stats::material_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_height_magnetization)          ::stats::height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);
               if(::stats::calculate_material_height_magnetization) ::stats::material_height_magnetization.calculate_magnetization(::atoms::x_spin_array, ::atoms::y_spin_array, ::atoms::z_spin_array, ::atoms::m_spin_array);

               // return before doing the GPU version
               return;
            }

            cu::stats::__update_stat (
                  cu::stats::system_mask,
                  cu::stats::system_magnetization,
                  cu::stats::system_mean_magnetization,
                  cu::stats::system_mask_size);

            cu::stats::__update_stat (
                  cu::stats::material_mask,
                  cu::stats::material_magnetization,
                  cu::stats::material_mean_magnetization,
                  cu::stats::material_mask_size);

            cu::stats::__update_stat (
                  cu::stats::height_mask,
                  cu::stats::height_magnetization,
                  cu::stats::height_mean_magnetization,
                  cu::stats::height_mask_size);

            cu::stats::__update_stat (
                  cu::stats::material_height_mask,
                  cu::stats::material_height_magnetization,
                  cu::stats::material_height_mean_magnetization,
                  cu::stats::material_height_mask_size);

            // increase the counter
            cu::stats::counter++;

         }

         void get ()
         {

            // If CPU stats calculation do nothing
            if(vcuda::internal::stats::use_cpu) return;

            cu::stats::__get_stat (
                  cu::stats::system_magnetization,
                  cu::stats::system_mean_magnetization,
                  ::stats::system_magnetization
                  );

            cu::stats::__get_stat (
                  cu::stats::material_magnetization,
                  cu::stats::material_mean_magnetization,
                  ::stats::material_magnetization
                  );

            cu::stats::__get_stat (
                  cu::stats::height_magnetization,
                  cu::stats::height_mean_magnetization,
                  ::stats::height_magnetization
                  );

            cu::stats::__get_stat (
                  cu::stats::material_height_magnetization,
                  cu::stats::material_height_mean_magnetization,
                  ::stats::material_height_magnetization
                  );
         }

         void reset ()
         {

            // reset magnetization statistics
            if(vcuda::internal::stats::use_cpu){
               if(::stats::calculate_system_magnetization)          ::stats::system_magnetization.reset_magnetization_averages();
               if(::stats::calculate_material_magnetization)        ::stats::material_magnetization.reset_magnetization_averages();
               if(::stats::calculate_height_magnetization)          ::stats::height_magnetization.reset_magnetization_averages();
               if(::stats::calculate_material_height_magnetization) ::stats::material_height_magnetization.reset_magnetization_averages();
               return;
            }

            cu::stats::counter = 0L;

            cu::stats::__reset_stat (
                  cu::stats::system_magnetization,
                  cu::stats::system_mean_magnetization
                  );

            cu::stats::__reset_stat (
                  cu::stats::material_magnetization,
                  cu::stats::material_mean_magnetization
                  );

            cu::stats::__reset_stat (
                  cu::stats::height_magnetization,
                  cu::stats::height_mean_magnetization
                  );

            cu::stats::__reset_stat (
                  cu::stats::material_height_magnetization,
                  cu::stats::material_height_mean_magnetization
                  );
         }

   } /* stats */

   namespace internal
   {
      namespace stats
      {
          typedef thrust::tuple<double, double> double2;

          // This functor implements the dot product between 3d vectors
          struct D2Product : public thrust::unary_function<double2, double>
          {
              __host__ __device__
                  double operator()(const double2& a) const
                  {
                      return thrust::get<0>(a) * thrust::get<1>(a);
                  }
          };

         void __update_stat (
               const cu_index_array_t & mask,
               cu_real_array_t & stat,
               cu_real_array_t & mean_stat,
               int mask_size
               )
         {

            if (mask_size < 1) return; // Nothing to do
            else if ( mask_size == 1)
            {
                typedef cu_real_array_t::iterator real_it;
                typedef thrust::tuple< real_it, real_it> real_it_tuple;
                typedef thrust::zip_iterator<real_it_tuple> real2_it;

                real2_it begin = thrust::make_zip_iterator( make_tuple( cu::atoms::x_spin_array.begin(), cu::atoms::spin_norm_array.begin() ) );
                real2_it end = thrust::make_zip_iterator( make_tuple( cu::atoms::x_spin_array.end(), cu::atoms::spin_norm_array.end() ) );
                stat[0] = thrust::transform_reduce( begin, end, D2Product(), 0.0, thrust::plus<float>()  );

                begin = thrust::make_zip_iterator( make_tuple( cu::atoms::y_spin_array.begin(), cu::atoms::spin_norm_array.begin() ) );
                end = thrust::make_zip_iterator( make_tuple( cu::atoms::y_spin_array.end(), cu::atoms::spin_norm_array.end() ) );
                stat[1] = thrust::transform_reduce( begin, end, D2Product(), 0.0, thrust::plus<float>()  );

                begin = thrust::make_zip_iterator( make_tuple( cu::atoms::z_spin_array.begin(), cu::atoms::spin_norm_array.begin() ) );
                end = thrust::make_zip_iterator( make_tuple( cu::atoms::z_spin_array.end(), cu::atoms::spin_norm_array.end() ) );
                stat[2] = thrust::transform_reduce( begin, end, D2Product(), 0.0, thrust::plus<float>()  );

                double ms = thrust::reduce( cu::atoms::spin_norm_array.begin(), cu::atoms::spin_norm_array.end() );
                double m = sqrt( stat[0]*stat[0] + stat[1]*stat[1] + stat[2]*stat[2]);
                stat[0] /= m;
                stat[1] /= m;
                stat[2] /= m;
                stat[3] = m / ms;

                mean_stat[0] += stat[0];
                mean_stat[1] += stat[1];
                mean_stat[2] += stat[2];
                mean_stat[3] += stat[3];

                return;
            }

            // Clean up the stat buffer
            thrust::fill(
                  stat.begin(),
                  stat.end(),
                  0.0);

            const int * d_mask = thrust::raw_pointer_cast (
                  mask.data());
            cu_real_t * d_stat = thrust::raw_pointer_cast (
                  stat.data());
            cu_real_t * d_accu = thrust::raw_pointer_cast (
                  mean_stat.data());

            cu_real_t * d_x_spin = thrust::raw_pointer_cast(
                  cu::atoms::x_spin_array.data());
            cu_real_t * d_y_spin = thrust::raw_pointer_cast(
                  cu::atoms::y_spin_array.data());
            cu_real_t * d_z_spin = thrust::raw_pointer_cast(
                  cu::atoms::z_spin_array.data());
            cu_real_t * d_spin_norm = thrust::raw_pointer_cast(
                  cu::atoms::spin_norm_array.data());

            int n_bins = mask_size;
            int n_atoms = mask.size ();

            if (n_bins < 8) {
               hist_by_key_smaller_mask <<< n_bins * 4, 512 >>> (
                     d_x_spin,
                     d_y_spin,
                     d_z_spin,
                     d_spin_norm,
                     d_mask,
                     d_stat,
                     n_bins,
                     n_atoms
                     );
               check_cuda_errors (__FILE__, __LINE__);
            }
            else if (n_bins < 128)
            {

                // Use the shared memory implementation

               int n_bytes = 4 * mask_size * sizeof(cu_real_array_t::value_type);
               hist_by_key_small_mask <<< cu::grid_size, cu::block_size, n_bytes >>> (
                     d_x_spin,
                     d_y_spin,
                     d_z_spin,
                     d_spin_norm,
                     d_mask,
                     d_stat,
                     n_bins,
                     n_atoms
                     );
               check_cuda_errors (__FILE__, __LINE__);
            }
            else
            {

               // Use the brute force implementation

               hist_by_key_big_mask <<< cu::grid_size, cu::block_size >>> (
                     d_x_spin,
                     d_y_spin,
                     d_z_spin,
                     d_spin_norm,
                     d_mask,
                     d_stat,
                     n_bins,
                     n_atoms
                     );
               check_cuda_errors (__FILE__, __LINE__);
            }

             // Reduce and accumulate

            int gs = n_bins / cu::block_size + 1;
            update_norm_and_accum <<< gs , cu::block_size >>> (
                  d_stat,
                  d_accu,
                  n_bins
                  );

            check_cuda_errors (__FILE__, __LINE__);

         }


         void __get_stat (
               const cu_real_array_t& stat,
               const cu_real_array_t& mean_stat,
               ::stats::magnetization_statistic_t& local_stat
               )
         {

            /*
             * Copy to local arrays
             */

            thrust::host_vector<cu_real_t> h_stat(stat.size());
            thrust::host_vector<cu_real_t> h_mean_stat(mean_stat.size());

            thrust::copy(stat.begin(), stat.end(), h_stat.begin());
            thrust::copy(mean_stat.begin(), mean_stat.end(), h_mean_stat.begin());

            /*
             * Call the method in the magnetization_statistic_t instance
             */

            std::vector<double> stl_stat (h_stat.begin(), h_stat.end());
            std::vector<double> stl_mean_stat (h_mean_stat.begin(), h_mean_stat.end());

            // We're keeping the averages in the GPU
            local_stat.reset_magnetization_averages();
            local_stat.set_magnetization (
                  stl_stat,
                  stl_mean_stat,
                  cu::stats::counter);
            check_cuda_errors (__FILE__, __LINE__);

         }


         void __reset_stat (
               cu_real_array_t& stat,
               cu_real_array_t& mean_stat
               )
         {
            thrust::fill(
                  stat.begin(),
                  stat.end(),
                  0.0);
            thrust::fill(
                  mean_stat.begin(),
                  mean_stat.end(),
                  0.0);
            check_cuda_errors (__FILE__, __LINE__);
         }

         __global__ void hist_by_key_smaller_mask (
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               )
         {
            cu_real_t val = 0; // Initialize the register

            int rol = blockIdx.x % 4;
            int bin = blockIdx.x / 4;
            // Block stride loop (each block reads the whole array)
            for ( int i = threadIdx.x;
                  i < n_atoms;
                  i += blockDim.x)
            {
               // Single spin test
               // if (i != 0) continue;
               cu_real_t mu_s = norm_spin[i];
               if (rol == 0) val += mask[i] == bin ? x_spin[i] * mu_s : 0;
               if (rol == 1) val += mask[i] == bin ? y_spin[i] * mu_s : 0;
               if (rol == 2) val += mask[i] == bin ? z_spin[i] * mu_s : 0;
               if (rol == 3) val += mask[i] == bin ? mu_s : 0;
            }

            typedef cub::BlockReduce<cu_real_t, 512> BlockReduce;
            __shared__ typename BlockReduce::TempStorage temp_storage;
            val = BlockReduce(temp_storage).Sum(val);

            if (threadIdx.x == 0) {
               hist[4 * bin + rol] += val;
            }
         }

         /**
           * Naivive implementation of histogram using shared memory
           * FIXME: Prone to improvements
           */
         __global__ void hist_by_key_small_mask (
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               )
         {
            extern __shared__ cu_real_t block_hist[];

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
               cu_real_t mu_s = norm_spin[i];
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
               const cu_real_t * __restrict__ x_spin,
               const cu_real_t * __restrict__ y_spin,
               const cu_real_t * __restrict__ z_spin,
               const cu_real_t * __restrict__ norm_spin,
               const int * __restrict__ mask,
               cu_real_t * hist,
               int n_bins,
               int n_atoms
               )
         {

            for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
                  i < n_atoms;
                  i += blockDim.x * gridDim.x)
            {

                // Store stuff in the main memory

               int bin = mask[i];
               cu_real_t mu_s = norm_spin[i];
               cu::atomicAdd (hist + 4 * bin + 0, x_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 1, y_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 2, z_spin[i] * mu_s);
               cu::atomicAdd (hist + 4 * bin + 3, mu_s);
            }
         }


         __global__ void update_norm_and_accum (
               cu_real_t * hist,
               cu_real_t * accum,
               int n_bins
               )
         {

            for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
                  i < n_bins;
                  i += blockDim.x * gridDim.x)
            {
               cu_real_t mx = hist[4 * i + 0];
               cu_real_t my = hist[4 * i + 1];
               cu_real_t mz = hist[4 * i + 2];
               cu_real_t ms = hist[4 * i + 3];

               if (ms > 0.0) {

                  // cuda overloads reciprocal sqrt for float and double
                  cu_real_t imm = rsqrt(mx * mx + my * my + mz * mz);

                  hist[4 * i + 0] = mx * imm;
                  hist[4 * i + 1] = my * imm;
                  hist[4 * i + 2] = mz * imm;
                  hist[4 * i + 3] = ms / imm;
                  accum[4 * i + 0] += mx * imm;
                  accum[4 * i + 1] += my * imm;
                  accum[4 * i + 2] += mz * imm;
                  accum[4 * i + 3] += ms / imm;
               } else {
                  // Just wipe the histogram
                  // FIXME Even this could be removed
                  for (int j = 0; j < 4; j++)
                     hist[4 * i + j] = 0.0;
               }
            }
         }

      } /* stats */
   } /* internal */

#endif

} // end of namespace cuda
