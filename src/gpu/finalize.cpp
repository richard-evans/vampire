//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2014. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"
#include "mgpu.hpp"
#include "errors.hpp"

namespace gpu{

   //----------------------------------------------------------------------------------
   // Function to call correct de-initialization function depending on the backend (cuda)
   //----------------------------------------------------------------------------------
   void finalize(){

      #ifdef CUDA_MULTIGPU
         mgpu::finalize();
      #endif

      #ifdef CUDA
         vcuda::finalize();
      #endif

      return;
   }

} // end of namespace gpu
