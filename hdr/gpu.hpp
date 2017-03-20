#ifndef GPU_H_
#define GPU_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015, 2016. All rights reserved.
//
//-----------------------------------------------------------------------------
#include <string>

namespace gpu{

   //-----------------------------------------------------------------------------
   // Variables used for GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool acceleration; // flag to enable gpu_acceleration
   extern bool cpu_stats; // flag to calculate stats using cpu

   extern int num_threads; // number of threads to run per kernel
   extern int platform; // which platform to use for acceleration (OpenCL)
   extern int device; // int specifying gpu device to use for simulation


   //-----------------------------------------------------------------------------
   // Functions for GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   extern void initialize();
   extern void llg_heun();
   extern void finalize();

   //-----------------------------------------------------------------------------
   // Functions for GPU configuration output
   //-----------------------------------------------------------------------------
   namespace config{
      extern void synchronise();
   }

   //-----------------------------------------------------------------------------
   // Functions for GPU statistics calculation
   //-----------------------------------------------------------------------------
   namespace stats{

      extern void update();
      extern void reset();
      extern void get();

   }

} // end of gpu namespace

#endif //GPU_H_
