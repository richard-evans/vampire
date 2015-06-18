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
      __global__ void init_rng (curandState * state, size_t seed)
      {
         size_t tid = blockIdx.x + blockDim.x + threadIdx.x;
         curand_init (seed, tid, 0, state + tid);
      }
   }

#endif
} /* vcuda */