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

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "internal.hpp"

namespace cu = cuda::internal;

namespace cuda {

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

#ifdef CUDA

      bool success = true;

      success = success || __initialize_atoms ();
      success = success || __initialize_fields ();
      success = success || __initialize_cells ();
      success = success || __initialize_materials ();
      success = success || __initialize_topology ();

      // send topology information

      // Successful initialization
      return success;

#endif

      // Default (initializtion failed)
      return false;

   }

} // end of namespace cuda
