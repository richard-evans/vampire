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

namespace cuda{

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

      #ifdef CUDA

         // Successful initialization
         return true;

      #endif

      // Default (initializtion failed)
      return false;

   }

} // end of namespace cuda
