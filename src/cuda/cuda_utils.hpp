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
#include <cub-1.5.2/cub/cub.cuh>

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

/**
 * This will butterfly reduce a value within a wrap, however
 * this won't work with doubles.
 */
template <typename T>
__inline__ __device__ T warpReduceSum (T val) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2)
        val += __shfl_down(val, offset);
    return val;
}

/**
 * This will warp reduce and shared memory reduce a value within a
 * block, this is pretty general but requires the warp reduce, these
 * essentials might be availabe in CUB.
 */
template <typename T>
__inline__ __device__ T blockReduceSum (T val) {
    static __shared__ T shared[32]; // Popamoly this is the max
    int laneId = threadIdx.x & 0x1f;
    int warpId = threadIdx.x / warpSize;
    val = warpReduceSum(val);
    if (laneId == 0) shared[warpId] = val;
    __syncthreads(); // Wait for al warp reductions
    // FIXME: This might need a + 1
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[laneId] : 0;
    if (warpId == 0) val = warpReduceSum(val);
    return val;
}

#endif
