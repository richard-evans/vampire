//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"
#include "errors.hpp"
#include "vio.hpp"
//#include "opencl.hpp"

namespace gpu{

   //-------------------------------------------------------------------------------
   // Function to call correct initialization function depending on cuda,opencl etc
   //-------------------------------------------------------------------------------
   void initialize(){

      // flag to check for successful initialization
      bool initialized=false;

       #ifdef CUDA
         initialized = cuda::initialize();
      #elseif OPENCL
         initialized = opencl::initialize();
      #endif

      // Check for no initialization
      if(initialized == false){
         std::cerr << "Error: no support for GPU acceleration, please recompile using -D CUDA compile" << std::endl;
         std::cerr << "time flag with cuda development libraries v6.0 or greater." << std::endl;
         std::cerr << "Aborting program." << std::endl;
         zlog << zTs() << "Error: no support for GPU acceleration, please recompile using -D CUDA compile time flag with cuda development libraries v6.0 or greater." << std::endl;
         err::vexit();
      }

      return;
   }

} // end of namespace gpu
