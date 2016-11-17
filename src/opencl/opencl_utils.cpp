#include "opencl_include.hpp"

#include <fstream>
#include <string>

#include "internal.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      // function to build OpenCL kernel from source file
      cl::Kernel build_kernel_from_file(const std::string &file_name,
                                        const std::string &kernel_name,
                                        const cl::Context &context,
                                        const cl::Device  &device,
                                        const std::string &opts="")
      {
         std::ifstream source_file(file_name);
         std::string source_code(std::istreambuf_iterator<char>(source_file),
                                (std::istreambuf_iterator<char>()));
         cl::Program::Sources source(1, std::make_pair(source_code.c_str(),
                                                       source_code.length()+1));


         cl::Program program = cl::Program(context, source);

         std::string prg_opts(opts);

#ifdef OPENCL_DP
         prg_opts.append(" -DOPENCL_DP");
#endif
#ifdef OPENCL_USE_NATIVE_FUNCTIONS
         prg_opts.append(" -DOPENCL_USE_NATIVE_FUNCTIONS");
#endif

         program.build({device}, prg_opts.c_str());

#ifdef OPENCL_DEBUG
         vcl::OCLLOG << "Building from " << file_name;
         vcl::OCLLOG << " using kernel " << kernel_name << std::endl;
         vcl::OCLLOG << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(device) << std::endl;
         vcl::OCLLOG << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
#endif

         cl::Kernel ret(program, kernel_name.c_str());
         return ret;
      }
   }
}

#endif // OPENCL
