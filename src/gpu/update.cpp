//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2021. All rights reserved.
//
//------------------------------------------------------------------------------
//
// C++ standard library headers

// Vampire headers
#include "gpu.hpp"
#include "cuda.hpp"

namespace gpu{

//------------------------------------------------------------------------------
// Wrapper function to update dipolar field with GPU accelleration
//------------------------------------------------------------------------------
void update_dipolar_fields(){

   #ifdef CUDA
      vcuda::update_dipolar_fields();
   #elif OPENCL
      // vopencl::update_dipolar_field(); TBD
   #endif

   return;

}
} // end of gpu namespace
