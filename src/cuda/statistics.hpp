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
      } /* stats */
   } /* internal */
   #endif
} /* vcuda */

#endif
