//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "opencl_include.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "errors.hpp"

#include "internal.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      void kernel_call(const cl::Kernel &k,            /* kernel to enqueue */
                       const cl::CommandQueue &q=vcl::queue,      /* into this queue */
                       const cl::NDRange gbl=vcl::global,          /* total number of work items */
                       const cl::NDRange lcl=vcl::local) noexcept /* number of work items in group */
      {
         const auto err = q.enqueueNDRangeKernel(k, cl::NullRange, gbl, lcl);
         if (err != CL_SUCCESS)
         {
            std::cerr << "Error enqueuing kernel " << k.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
            std::cerr << "Error code " << err << std::endl;
            ::err::vexit();
         }
      }

      // function to build OpenCL kernel from source file
      cl::Kernel build_kernel_from_file(const std::string &file_name,
                                        const std::string &kernel_name,
                                        const std::string &opts="",
                                        const cl::Device  &device=vcl::default_device,
                                        const cl::Context &context=vcl::context) noexcept
      {
         std::ifstream source_file(file_name);
         if (!source_file)
         {
            std::cerr << "File " << file_name << " not found!" << std::endl;
         }
         std::string source_code(std::istreambuf_iterator<char>(source_file),
                                 (std::istreambuf_iterator<char>()));
         cl::Program::Sources source(1, std::make_pair(source_code.c_str(),
                                                       source_code.length()+1));

         cl::Program program = cl::Program(context, source);

         cl_int err = program.build({device}, opts.c_str());

#ifdef OPENCL_LOG
         vcl::OCLLOG << "Building from " << file_name;
         vcl::OCLLOG << " using kernel " << kernel_name << std::endl;
         vcl::OCLLOG << "using options ";
         vcl::OCLLOG << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device) << std::endl;
         vcl::OCLLOG << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
#endif

         if (err != CL_SUCCESS)
         {
            std::cerr << "Error building kernel " << kernel_name << " in file " << file_name << std::endl;
            std::cerr << "cl::Program::build() returned " << err << std::endl;
            ::err::vexit();
         }

         cl::Kernel ret(program, kernel_name.c_str());

         return ret;
      }
   }
}

#endif // OPENCL
