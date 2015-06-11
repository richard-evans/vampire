//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

#include <vector>

// Thrust library headers

#include <thrust/device_vector.hpp>
#include <thrust/copy.hpp>

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "internal.hpp"

namespace cuda{

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

#ifdef CUDA

      // send the spin information
      // send the material information
      // send the macro-cell information

      // Successful initialization
      return true;

#endif

      // Default (initializtion failed)
      return false;

   }

} // end of namespace cuda
