// defines constants, types and functions to be used in
// OpenCL source files

#ifndef VOPENCL_CL_DEFS_H_
#define VOPENCL_CL_DEFS_H_

#define PI 3.14159265358979323846

#ifdef OPENCL_DP

// double precision (64-bit) real and unsigned integer types
typedef double real_t;
typedef ulong  uint_t;

// double precision (64-bit) function
#define ATOMIC_CMPXCHG(a,b,c) atom_cmpxchg(a,b,c)

#else

// single precision (32-bit) real and unsigned interge types
typedef float real_t;
typedef uint  uint_t;

// single precision (32-bit) function
#define ATOMIC_CMPXCHG(a,b,c) atomic_cmpxchg(a,b,c)

#endif // OPENCL_DP



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
