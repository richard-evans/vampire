//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_INTERNAL_H_
#define VOPENCL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the vopencl module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <chrono>
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

      namespace time
      {
         typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_t;
         extern time_t sim_start;

#ifdef OPENCL_TIME_KERNELS
         extern double spin_fields;
         extern double mat_mul;
         extern double rng;
         extern double external_fields;
         extern double predictor_step;
         extern double corrector_step;
#endif // OPENCL_TIME_KERNELS
      }

      // number of work items
      extern cl::NDRange global;
      extern cl::NDRange global_other;

      // let implementation choose work group size
      const cl::NDRange local(cl::NullRange);

      extern cl::Context context;
      extern cl::Context context_other;
      extern cl::Device default_device;
      extern cl::Device extra_device;
      extern cl::CommandQueue queue;
      extern cl::CommandQueue queue_other;

#ifdef OPENCL_LOG
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
