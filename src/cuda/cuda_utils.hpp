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


//---------------------------------------------------------------------
// Defines cuda utility functions
// Author: Matt Ellis
// Created: 18/6/2015
//---------------------------------------------------------------------


#ifndef __CUDA_UTILS_HPP__
#define __CUDA_UTILS_HPP__


inline void check_cuda_errors(const char *filename, const int line_number)
{
#ifdef CODE_DEBUG
    cudaThreadSynchronize();
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
        exit(-1);
    }
#endif
}



#endif
