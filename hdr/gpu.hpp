#ifndef GPU_H_
#define GPU_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

namespace gpu{

   //-----------------------------------------------------------------------------
   // Variables used for GPU acceleration
   //-----------------------------------------------------------------------------
   extern bool acceleration; // flag to enable gpu_acceleration

   //-----------------------------------------------------------------------------
   // Functions for GPU acceleration
   //-----------------------------------------------------------------------------
   extern void initialize();
   extern void llg_heun();
   extern void stats_update();
   extern void finalize();

} // end of gpu namespace

#endif //GPU_H_
