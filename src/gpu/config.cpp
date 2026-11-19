//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"

namespace gpu{

   namespace config{

      //-------------------------------------------------------------------------------
      // Function to synchronise spin array from device to cpu
      //-------------------------------------------------------------------------------
      void synchronise(){

         // check for gpu functions and data initialised
         if( !gpu::initialized ) return;

         #ifdef CUDA
            vcuda::config::synchronise();
         #endif

         return;
      }

   } // end of config namespace

} // end of namespace gpu
