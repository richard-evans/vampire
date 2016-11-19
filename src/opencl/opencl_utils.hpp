#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <fstream>
#include <string>

#include "opencl_include.hpp"
#include "internal.hpp"

namespace vopencl
{
   namespace internal
   {
      cl::Kernel build_kernel_from_file(const std::string &file_name,
                                        const std::string &kernel_name,
                                        const cl::Context &context,
                                        const cl::Device &device,
                                        const std::string &opts="");

      // class used to execute kernels
      class kernel_call
      {
         cl::Kernel clk;
         unsigned i;

         void pass_args(void){}

         template <typename Car, typename... Cdr>
         void pass_args(Car car, Cdr... cdr)
         {
            clk.setArg(i++, car);
            this->pass_args(cdr...);
         }

      public:

         template <typename... Ts>
         kernel_call(const cl::Kernel &k,
                     const cl::CommandQueue &q,
                     const cl::NDRange gbl,
                     const cl::NDRange lcl,
                     Ts... Args) :
            clk(k), i(0)
         {
            this->pass_args(Args...);
            q.enqueueNDRangeKernel(clk, cl::NullRange, gbl, lcl);
         }
      };
   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
