#ifndef CUDA_H_
#define CUDA_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

namespace cuda{

#ifdef CUDA

   //-----------------------------------------------------------------------------
   // Variables used for cuda GPU acceleration
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Functions for cuda GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool initialize();
   extern void llg_heun();
   extern void stats_update();
   extern void finalize();

#endif

} // end of gpu namespace

#endif //CUDA_H_
