#ifndef VOPENCL_OPENCL_INCLUDE_HPP_
#define VOPENCL_OPENCL_INCLUDE_HPP_

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif

// if debugging then buffers which are normally not read from the host
// will be so need to redefine these flags, otherwise undefined behaviour
#ifdef OPENCL_DEBUG
#undef CL_MEM_HOST_WRITE_ONLY
#define CL_MEM_HOST_WRITE_ONLY 0

#undef CL_MEM_HOST_NO_ACCESS
#define CL_MEM_HOST_NO_ACCESS 0
#endif // OPENCL_DEBUG

#endif // VOPENCL_OPENCL_INCLUDE_HPP_
