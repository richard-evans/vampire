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
#include "internal.hpp"

namespace cuda{

   //--------------------------------------------------------------------------
   // Function to perform a single heun step
   //--------------------------------------------------------------------------
   void llg_heun(){

      #ifdef CUDA
      /* set up and call the kernels */
      /* assume that you have the data already
       * in the device */

      llg_heun_first_kernel <<< thin, thang >>> (
            /*
               cuda::internal::spins,
               cuda::internal::other_stuff,
               ...
            */
            );

      #endif

      return;
   }

} // end of namespace cuda
