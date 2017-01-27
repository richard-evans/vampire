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

      // number of work items
      extern cl::NDRange global;

      // let implementation choose work group size
      const cl::NDRange local(cl::NullRange);

      extern cl::Context context;
      extern cl::Device default_device;
      extern cl::CommandQueue queue;

#ifdef OPENCL_DEBUG
      extern std::ofstream OCLLOG;
#endif
      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      bool initialize_atoms(void) noexcept;
      bool initialize_fields(void) noexcept;
      bool initialize_cells(void) noexcept;
      bool initialize_materials(void) noexcept;
      bool initialize_topology(void) noexcept;
      bool initialize_stats(void) noexcept;
      bool initialize_rng(void) noexcept;

      void update_external_fields() noexcept;
      void update_dipolar_fields() noexcept;
      void update_cell_magnetizations() noexcept;

      void finalize(void) noexcept;
      
   } // end of internal namespace

} // end of vopencl namespace

#endif //VOPENCL_INTERNAL_H_
