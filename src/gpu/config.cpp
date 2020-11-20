//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"
#include "vopencl.hpp"

namespace gpu{

   namespace config{

      //-------------------------------------------------------------------------------
      // Function to synchronise spin array from device to cpu
      //-------------------------------------------------------------------------------
      void synchronise(){

         #ifdef CUDA
            vcuda::config::synchronise();
         #elif OPENCL
            vopencl::config::synchronise();
         #endif

         return;
      }

   } // end of config namespace

} // end of namespace gpu
