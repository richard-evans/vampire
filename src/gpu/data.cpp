//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "gpu.hpp"

namespace gpu{

   //-----------------------------------------------------------------------------
   // Shared variables used for GPU acceleration
   //-----------------------------------------------------------------------------
   // default is on for GPU compiled code
   #ifdef CUDA
      bool acceleration = true; // flag to enable gpu_acceleration
   #elif OPENCL
      bool acceleration = true; // flag to enable gpu_acceleration
   #else
      bool acceleration = false; // flag to enable gpu_acceleration
   #endif

   bool cpu_stats = false; // flag to calculate stats using cpu

   int num_threads = 0; // number of threads to use per kernel

   //-----------------------------------------------------------------------------
   // Shared data structures for statistics calculation
   //-----------------------------------------------------------------------------
   namespace stats{

      long counter; // saves the number of averages accumulated on the gpu

      // temporary arrays for storing system magnetization data from the gpu
      std::vector<double> system_magnetization(0);
      std::vector<double> system_mean_magnetization(0);

      // temporary arrays for storing material magnetization data from the gpu
      std::vector<double> material_magnetization(0);
      std::vector<double> material_mean_magnetization(0);

      // temporary arrays for storing height magnetization data from the gpu
      std::vector<double> height_magnetization(0);
      std::vector<double> height_mean_magnetization(0);

      // temporary arrays for storing height-material magnetization data from the gpu
      std::vector<double> material_height_magnetization(0);
      std::vector<double> material_height_mean_magnetization(0);

   } // end of stats namespace

   namespace internal{

   } // end of internal namespace

} // end of gpu namespace
