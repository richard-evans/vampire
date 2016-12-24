#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <fstream>
#include <iostream>
#include <string>

#include "errors.hpp"

#include "internal.hpp"
#include "opencl_include.hpp"

#define UNUSED(x) (void)(x)

static std::string get_error(cl_int err)
{
   switch (err)
   {
   case CL_INVALID_VALUE:
      return "CL_INVALID_VALUE";
   case CL_INVALID_KERNEL:
      return "CL_INVALID_KERNEL";
   case CL_OUT_OF_RESOURCES:
      return "CL_OUT_OF_RESOURCES";
   case CL_OUT_OF_HOST_MEMORY:
      return "CL_OUT_OF_HOST_MEMORY";
   default:
      return "unknown error code";
   }
}

// functions to pass given arguments to kernel k which recursively sets argument i
// the following function is needed for when all the arguments have been set
static void pass_args(cl::Kernel &k, unsigned i)
{
   UNUSED(k);
   UNUSED(i);
}

// function which takes an arg from args and gives it to the kernel for arg i
template <typename Car, typename... Cdr>
static void pass_args(cl::Kernel &k, unsigned i, Car car, Cdr... cdr)
{
   cl_int err;
   if ((err = k.setArg(i, car)) == CL_SUCCESS)
      pass_args(k, i+1, cdr...);
   else
   {
      std::cerr << "Error setting kernel argument " << (i-1);
      std::cerr << " in kernel " << k.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
      std::cerr << "error code " << get_error(err) << std::endl;
      ::err::vexit();
   }
}

namespace vopencl
{
   namespace internal
   {
      // function to return a kernel from a given file name and kernel name
      // the opts string can be used to define things as with gcc
      // for example "-DN" defines N (for use in #ifdef)
      // and "-DN=10" will replace N with 10 in the kernel code
      cl::Kernel build_kernel_from_file(const std::string &file_name,
                                        const std::string &kernel_name,
                                        const cl::Context &context,
                                        const cl::Device  &device,
                                        const std::string &opts="");

      // function to enqueue a kernel
      // the OpenCL C++ spec for this seems pretty fluid so it's probably better not to use it
      template <typename... Ts>
      static void kernel_call(cl::Kernel &k,        /* kernel to enqueue */
                              cl::CommandQueue &q,  /* into this queue */
                              cl::NDRange gbl,      /* total number of work items */
                              cl::NDRange lcl,      /* number of work items in group */
                              Ts... Args)           /* kernel arguments */
      {
         pass_args(k, 0, Args...);
         auto err = q.enqueueNDRangeKernel(k, cl::NullRange, gbl, lcl);
         if (err != CL_SUCCESS)
         {
            std::cerr << "Error enqueuing kernel " << k.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
            std::cerr << "Error code " << err << std::endl;
            ::err::vexit();
         }
      }
   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
