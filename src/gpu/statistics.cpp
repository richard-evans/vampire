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

   namespace stats{

      //-------------------------------------------------------------------------------
      // Function to call correct statistics update function on device
      //-------------------------------------------------------------------------------
      void update(){

         #ifdef CUDA
            vcuda::stats::update();
         #elif OPENCL
            vopencl::stats::update();
         #endif

         return;
      }

      //-------------------------------------------------------------------------------
      // Function to get statistics module data from device for output
      //
      // Remember to use the set_magnetization in the magnetization_statistic_t instances
      //
      // There is a performance penalty for rapid updates of statistics for output
      // to disk.
      //
      //-------------------------------------------------------------------------------
      void get(){

         #ifdef CUDA
            vcuda::stats::get();
         #elif OPENCL
            vopencl::stats::get();
         #endif

         return;
      }

      //-------------------------------------------------------------------------------
      // Function to reset statistics counters on device
      //-------------------------------------------------------------------------------
      void reset(){

         #ifdef CUDA
            vcuda::stats::reset();
         #elif OPENCL
            vopencl::stats::reset();
         #endif

         return;
      }

   } // end of stats namespace

} // end of namespace gpu
