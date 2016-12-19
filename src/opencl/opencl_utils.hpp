#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <fstream>
#include <iostream>
#include <string>

#include "errors.hpp"

#include "internal.hpp"
#include "opencl_include.hpp"

#define UNUSED(x) (void)(x)

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
   if ((err = k.setArg(i++, car)) == CL_SUCCESS)
      pass_args(k, i, cdr...);
   else
   {
      std::cerr << "Error setting kernel argument " << (i-1);
      std::cerr << " in kernel " << k.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
      std::cerr << "error code " << err << std::endl;
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
         q.enqueueNDRangeKernel(k, cl::NullRange, gbl, lcl);
      }
   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
