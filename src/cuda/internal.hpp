#ifndef CUDA_INTERNAL_H_
#define CUDA_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// cuda implementation. These functions should 
// not be accessed outside of the local temperature pulse code.
//---------------------------------------------------------------------
namespace cuda{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the cuda implementation
      //-----------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // Shared functions and kernels used for the cuda implementation
      //-----------------------------------------------------------------------------

      __global__ void llg_heun_first_kernel (/* whole bunch of arrays */);

   } // end of iternal namespace
} // end of cuda namespace

#endif //CUDA_INTERNAL_H_
