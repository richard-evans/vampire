/**
 * @brief this file provides definitions for the off-topic internal function
 *        definitions.
 */

#include "data.hpp"
#include "internal.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda
{
#ifdef CUDA

   namespace internal
   {

      int block_size(256UL);
      int grid_size(512UL);

      __device__ double atomicAdd (double * address, double value)
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

      __global__ void init_rng (curandState * state, int seed)
      {
         int tid = blockIdx.x + blockDim.x + threadIdx.x;
         curand_init (seed, tid, 0, state + tid);
      }
   }

#endif
} /* vcuda */
