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
      int grid_size(32UL);

   }

#endif
} /* vcuda */
