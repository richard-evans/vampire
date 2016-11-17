//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------
//

#ifndef VOPENCL_INTERNAL_H_
#define VOPENCL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the vopencl module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>
#include <fstream>

// Vampire headers
#include "opencl_include.hpp"

namespace vopencl
{

   namespace internal
   {

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern cl::Context context;
      extern cl::Device default_device;
      
#ifdef OPENCL_DEBUG
      extern std::ofstream OCLLOG;
#endif
      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      bool initialize_atoms(void);
      bool initialize_fields(void);
      bool initialize_cells(void);
      bool initialize_materials(void);
      bool initialize_topology(void);
      bool initialize_stats(void);
      bool initialize_rng(void);

      void finalize(void);
      
   } // end of internal namespace

} // end of vopencl namespace

#endif //VOPENCL_INTERNAL_H_
