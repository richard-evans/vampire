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

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda {

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

#ifdef CUDA

      bool success = true;

      success = success || cu::__initialize_atoms ();
      success = success || cu::__initialize_fields ();
      success = success || cu::__initialize_cells ();
      success = success || cu::__initialize_materials ();
      success = success || cu::__initialize_topology ();

      // Successful initialization
      return success;

#else
      // Default (initializtion failed)
      return false;
#endif
   }

   bool finalize()
   {
#ifdef CUDA
      return cu::__finalize();
#else
      return false;
#endif
   }

} // end of namespace cuda
