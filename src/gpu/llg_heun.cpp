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
#include "errors.hpp"

namespace gpu{

   //--------------------------------------------------------------------------
   // Function to call correct llg_heun function depending on the backend (cuda)
   //--------------------------------------------------------------------------
   void llg_heun(){

      #ifdef CUDA
         vcuda::llg_heun();
      #endif

      return;
   }

} // end of namespace gpu
