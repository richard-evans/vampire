#ifndef CUDA_STATISTICS_H_
#define CUDA_STATISTICS_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) O Arbelaez Echeverri 2015. All rights reserved.
//
//------------------------------------------------------------------------------

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
