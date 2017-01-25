#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "errors.hpp"

#include "internal.hpp"
#include "opencl_include.hpp"

#define UNUSED(x) (void)(x)

static std::string get_error(const cl_int err) noexcept
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
      std::ostringstream code;
      code << err;
      return std::string("unknown error code ").append(code.str());
   }
}

// functions to pass given arguments to kernel k which recursively sets argument i
// the following function is needed for when all the arguments have been set
static void pass_args(const cl::Kernel &k, const unsigned i) noexcept
{
   UNUSED(k);
   UNUSED(i);
}

// function which takes an arg from args and gives it to the kernel for arg i
template <typename Car, typename... Cdr>
static void pass_args(cl::Kernel &k, const unsigned i, const Car &car, const Cdr &... cdr) noexcept
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
                                        const std::string &opts="") noexcept;

      template <typename... Ts>
      static void set_kernel_args(cl::Kernel &k, Ts &&... Args) noexcept
      {
          pass_args(k, 0, Args...);
      }

      // function to enqueue a kernel
      // the OpenCL C++ spec for this seems pretty fluid so it's probably better not to use it
      void kernel_call(const cl::Kernel &k,              /* kernel to enqueue */
                       const cl::CommandQueue &q,  /* into this queue */
                       const cl::NDRange gbl,      /* total number of work items */
                       const cl::NDRange lcl) noexcept;     /* number of work items in group */

   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
