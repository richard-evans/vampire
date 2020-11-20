//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// defines constants, types and functions to be used in
// OpenCL source files

#ifndef VOPENCL_CL_DEFS_H_
#define VOPENCL_CL_DEFS_H_

#define PI 3.14159265358979323846


//// floating point precision
#ifdef OPENCL_DP

// double precision (64-bit) real and unsigned integer types
typedef double  real_t;
typedef double2 real_t2;
typedef double3 real_t3;
typedef double4 real_t4;
typedef ulong   uint_t;

// double precision (64-bit) function
#define ATOMIC_CMPXCHG(a,b,c) atom_cmpxchg(a,b,c)
#define NORM(vec) normalize(vec)

#else // use single precision

// single precision (32-bit) real and unsigned integer types
typedef float  real_t;
typedef float2 real_t2;
typedef float3 real_t3;
typedef float4 real_t4;
typedef uint   uint_t;

// single precision (32-bit) function
#define ATOMIC_CMPXCHG(a,b,c) atomic_cmpxchg(a,b,c)
#define NORM(vec) fast_normalize(vec)

#endif // OPENCL_DP


//// data storage type
#ifdef OPENCL_USE_VECTOR_TYPE

// type used to access certain data
typedef real_t3 vec_t;

// each element of the array is real_t3 so can access as normal
#define VEC_LOAD(array, idx) array[idx]
#define VEC_STORE(array, idx, value) array[idx] = value
#define VEC_INCR(array, idx, value) array[idx] += value

#else

// arrays stored in x,y,z,x,y,z format using ordinary floats/doubles
typedef real_t vec_t;

// use intrinsic functions to access each set of elements
#define VEC_LOAD(array, idx) vload3(idx, array)
#define VEC_STORE(array, idx, value) vstore3(value, idx, array)
#define VEC_INCR(array, idx, value) vstore3(vload3(idx, array)+value, idx, array)

#endif // OPENCL_USE_VECTOR_TYPE


// native functions may be significantly faster
// but may be less accurate
#ifdef OPENCL_USE_NATIVE_FUNCTIONS

#define POW(x,y) native_powr(x, y)
#define RSQRT(x) native_rsqrt(x)
#define  SQRT(x) native_sqrt(x)
#define   LOG(x) native_log(x)
#define   SIN(x) native_sin(x)
#define   COS(x) native_cos(x)

#else

#define POW(x,y) pow(x,y)
#define RSQRT(x) rsqrt(x)
#define  SQRT(x) sqrt(x)
#define   LOG(x) log(x)
#define   SIN(x) sin(x)
#define   COS(x) cos(x)

#endif // OPENCL_USE_NATIVE_FUNCTIONS

#endif // VOPENCL_CL_DEFS_H_
