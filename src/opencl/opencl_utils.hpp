//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_OPENCL_UTILS_HPP_
#define VOPENCL_OPENCL_UTILS_HPP_

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "errors.hpp"

#include "debug.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"

#define UNUSED(x) (void)(x)

namespace vcl = ::vopencl::internal;

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
   cl_int ec;
   if ((ec = k.setArg(i, car)) == CL_SUCCESS)
      pass_args(k, i+1, cdr...);
   else
   {
      std::cerr << "Error setting kernel argument " << (i-1);
      std::cerr << " in kernel " << k.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
      std::cerr << "error code " << ec << std::endl;
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
                                        const std::string &opts,
                                        const cl::Device  &device=vcl::default_device,
                                        const cl::Context &context=vcl::context) noexcept;


      template <typename... Ts>
      static void set_kernel_args(cl::Kernel &k, Ts &&... Args) noexcept
      {
         pass_args(k, 0, Args...);
      }

      // function to enqueue a kernel
      // the OpenCL C++ spec for this seems pretty fluid so it's probably better not to use it
      void kernel_call(const cl::Kernel &k,              /* kernel to enqueue */
                       const cl::CommandQueue &q=vcl::queue,  /* into this queue */
                       const cl::NDRange gbl=vcl::global,      /* total number of work items */
                       const cl::NDRange lcl=vcl::local) noexcept;     /* number of work items in group */

      // function to create a device buffer from a host vector
      // defaults to a read/write buffer written with a non blocking write
      template <typename T>
      cl::Buffer create_device_buffer(const std::vector<T> &host_vector,
                                      const cl_mem_flags mem_flags=CL_MEM_READ_WRITE,
                                      const cl_bool blocking=CL_FALSE,
                                      const cl::CommandQueue &queue=vcl::queue,
                                      const cl::Context &context=vcl::context)
      {
         const size_t buffer_size = host_vector.size() * sizeof(T);
         cl::Buffer device_buffer(context, mem_flags, buffer_size);

         if (buffer_size != 0)
         {
            cl_int ec = queue.enqueueWriteBuffer(device_buffer, blocking, 0, buffer_size, host_vector.data());

            if (ec != CL_SUCCESS)
            {
               std::cerr << "Error in create_device_buffer, enqueueWriteBuffer returned " << ec << std::endl;
               ::err::vexit();
            }

#ifdef OPENCL_DEBUG
            queue.finish();
            vcl::debug::verify_copy(device_buffer, host_vector);
#endif // OPENCL_DEBUG
         }
         return device_buffer;
      }
   }
}

#endif // VOPENCL_OPENCL_UTILS_HPP_
