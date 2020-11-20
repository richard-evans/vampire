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

      int block_size = 128;
      int grid_size = 32;
      timer cuda_timer;
   }

#endif
} /* vcuda */
