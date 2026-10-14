//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) O Arbelaez Echeverri 2015. All rights reserved.
//
//------------------------------------------------------------------------------

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
