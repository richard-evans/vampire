#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <fstream>
#include <iostream>
#include <string>

#include "errors.hpp"

#include "internal.hpp"
#include "opencl_include.hpp"

#define UNUSED(x) (void)(x)

static void pass_args(cl::Kernel &k, unsigned i)
{
   UNUSED(k);
   UNUSED(i);
}

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
      cl::Kernel build_kernel_from_file(const std::string &file_name,
                                        const std::string &kernel_name,
                                        const cl::Context &context,
                                        const cl::Device &device,
                                        const std::string &opts="");

      template <typename... Ts>
      static void kernel_call(cl::Kernel &k,
                              cl::CommandQueue &q,
                              cl::NDRange gbl,
                              cl::NDRange lcl,
                              Ts... Args)
      {
         pass_args(k, 0, Args...);
         q.enqueueNDRangeKernel(k, cl::NullRange, gbl, lcl);
      }
   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
