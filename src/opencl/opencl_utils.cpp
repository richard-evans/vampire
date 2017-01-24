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
                       const cl::CommandQueue &q,      /* into this queue */
                       const cl::NDRange gbl,          /* total number of work items */
                       const cl::NDRange lcl) noexcept /* number of work items in group */
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
                                        const cl::Context &context,
                                        const cl::Device  &device,
                                        const std::string &opts="") noexcept
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

         std::string prg_opts(opts);

         prg_opts.append(" -Isrc/opencl/cl"
#ifdef OPENCL_DP
                         " -DOPENCL_DP"
#endif
#ifdef OPENCL_USE_NATIVE_FUNCTIONS
                         " -DOPENCL_USE_NATIVE_FUNCTIONS"
#endif
#ifdef USE_VECTOR_TYPE
                         " -DUSE_VECTOR_TYPE"
#endif
            );


         cl_int err = program.build({device}, prg_opts.c_str());

#ifdef OPENCL_DEBUG
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
