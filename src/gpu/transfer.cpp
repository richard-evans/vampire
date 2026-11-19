//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "gpu.hpp"
#include "cuda.hpp"

namespace gpu{

//------------------------------------------------------------------------------
// Wrapper function to transfer spin data from GPU to CPU
//------------------------------------------------------------------------------
void transfer_spin_positions_from_gpu_to_cpu(){

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::transfer_spin_positions_from_gpu_to_cpu();
   #endif

   return;

}

//------------------------------------------------------------------------------
// Wrapper function to invalidate spin data on the CPU when it is not kept up
// to date by the GPU (statistics calculated on the GPU)
//------------------------------------------------------------------------------
void invalidate_host_spins(){

   if(err::check==true) std::cout << "gpu::invalidate_host_spins has been called" << std::endl;

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::invalidate_host_spins();
   #endif

   return;

}

//------------------------------------------------------------------------------
// Wrapper function to copy the spins to the host after integration when the
// program may read them there directly
//------------------------------------------------------------------------------
void update_host_spins(){

   if(err::check==true) std::cout << "gpu::update_host_spins has been called" << std::endl;

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::update_host_spins();
   #endif

   return;

}

//------------------------------------------------------------------------------
// Wrapper function to transfer dipole field data from CPU to GPU
//------------------------------------------------------------------------------
void transfer_dipole_fields_from_cpu_to_gpu(){

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::transfer_dipole_fields_from_cpu_to_gpu();
   #endif

   return;

}

//------------------------------------------------------------------------------
// Wrapper function to transfer dipole field data from GPU to CPU
//------------------------------------------------------------------------------
void transfer_dipole_cells_fields_from_gpu_to_cpu(){

   // check for gpu functions and data initialised
   if( !gpu::initialized ) return;

   #ifdef CUDA
      vcuda::transfer_dipole_cells_fields_from_gpu_to_cpu();
   #endif

   return;

}

} // end of gpu namespace
