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
#include <iostream>

namespace gpu{

//------------------------------------------------------------------------------
// Wrapper function to update dipolar field with GPU accelleration
//------------------------------------------------------------------------------
void update_tensor_dipolar_fields(){

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::update_dipolar_fields();
      std::cout << "update fields" <<std::endl;

   #elif OPENCL
      // vopencl::update_dipolar_field(); TBD
   #endif

   return;

}
void update_hierarchical_dipolar_fields(){
  std::cout << "update_hierarchical_dipolar_fields" << "\t" << gpu::initialized<< std::endl;
   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;
  std::cout << "update_hierarchical_dipolar_fields" <<std::endl;
   #ifdef CUDA
      vcuda::update_hierarchical_dipolar_fields();
      std::cout << "update fields" <<std::endl;
   #elif OPENCL
      // vopencl::update_dipolar_field(); TBD
   #endif

   return;

}
} // end of gpu namespace
