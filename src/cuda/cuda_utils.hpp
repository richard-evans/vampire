//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//
//
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Defines cuda utility functions
// Author: Matt Ellis
// Created: 18/6/2015
//-----------------------------------------------------------------------------


#ifndef __CUDA_UTILS_HPP__
#define __CUDA_UTILS_HPP__

#include <stdio.h>

//-----------------------------------------------------------------------------
// Wrapper function to test for errors in kernels if the debug flag is set
// Usage: check_cuda_errors( __FILE__, __LINE__)
// __FILE__ and __LINE__ defaults macros will input the position of the error
//-----------------------------------------------------------------------------

inline void check_cuda_errors(const char *filename, const int line_number)
{
#ifdef CUDA_DEBUG
   cudaThreadSynchronize();
   cudaError_t error = cudaGetLastError();
   if(error != cudaSuccess)
   {
      printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
      exit(-1);
   }
#endif
}

inline void check_device_memory( const char* filename, const int line_number)
{
#ifdef CUDA_DEBUG
   size_t avail;
   size_t total;
   cudaMemGetInfo( &avail, &total);
   size_t used = total - avail;
   printf( "CUDA device memory usage at %s:%i: Used: %f Mb, Free %f Mb\n", filename, line_number, float(used)/(1024*1024), float(avail)/(1024*1024));
#endif
}

#endif
