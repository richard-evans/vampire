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

// Vampire headers
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
   #elif OPENCL
      // vopencl::transfer_spin_positions_from_gpu_to_cpu(); TBD
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
   #elif OPENCL
      // vopencl::transfer_dipole_fields_from_cpu_to_gpu(); TBD
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
   #elif OPENCL
      // vopencl::transfer_dipole_cells_fields_from_gpu_to_cpu(); TBD
   #endif

   return;

}

} // end of gpu namespace
