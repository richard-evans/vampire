//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"
#include "errors.hpp"
#include "vopencl.hpp"

namespace gpu{

   //----------------------------------------------------------------------------------
   // Function to call correct de-initialization function depending on cuda,opencl etc
   //----------------------------------------------------------------------------------
   void finalize(){

      #ifdef CUDA
         vcuda::finalize();
      #elif OPENCL
         vopencl::finalize();
      #endif

      return;
   }

} // end of namespace gpu
