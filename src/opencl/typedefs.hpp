#ifndef VOPENCL_TYPEDEFS_HPP_
#define VOPENCL_TYPEDEFS_HPP_

#include <vector>

#include "internal.hpp"
#include "opencl_include.hpp"
#include "opencl_utils.hpp"

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
#ifdef OPENCL_DP
      typedef cl_double vcl_real_t;
#else
      typedef cl_float vcl_real_t;
#endif // OPENCL_DP


      // Class for set of three buffers, i.e. for x,y,z arrays
      // uses non blocking reads/writes
      class Buffer3D
      {
         std::vector<cl::Buffer> buffers;
         size_t buffer_size;

      public:

         Buffer3D(void) noexcept : buffers(3), buffer_size(0) {}

         // initialize without writing, but with size
         // e.g. for use when generating buffers on device
         Buffer3D(const cl::Context &c,
                  const cl_mem_flags fs,
                  const size_t size) noexcept
            : buffers(3), buffer_size(size)
         {
            buffers[0] = cl::Buffer(c, fs, size);
            buffers[1] = cl::Buffer(c, fs, size);
            buffers[2] = cl::Buffer(c, fs, size);
         }

          // initialize with 3 vectors to write data to device
         template <typename T>
         Buffer3D(const cl::Context &c,
                  const cl::CommandQueue &q,
                  const cl_mem_flags fs,
                  const std::vector<T> &xs,
                  const std::vector<T> &ys,
                  const std::vector<T> &zs) noexcept
            : buffers(3)
         {
            buffer_size = xs.size();

            buffers[0] = cl::Buffer(c, fs, sizeof(T)*xs.size());
            buffers[1] = cl::Buffer(c, fs, sizeof(T)*ys.size());
            buffers[2] = cl::Buffer(c, fs, sizeof(T)*zs.size());

            q.enqueueWriteBuffer(buffers[0], CL_FALSE, 0, sizeof(T)*xs.size(), &xs[0]);
            q.enqueueWriteBuffer(buffers[1], CL_FALSE, 0, sizeof(T)*ys.size(), &ys[0]);
            q.enqueueWriteBuffer(buffers[2], CL_FALSE, 0, sizeof(T)*zs.size(), &zs[0]);
         }

         // reads data from device, assumes host vectors already have enough capacity
         template <typename T>
         void copy_to_host(const cl::CommandQueue &q,
                           std::vector<T> &xs,
                           std::vector<T> &ys,
                           std::vector<T> &zs) const noexcept
         {
            q.enqueueReadBuffer(buffers[0], CL_FALSE, 0, sizeof(T)*xs.size(), &xs[0]);
            q.enqueueReadBuffer(buffers[1], CL_FALSE, 0, sizeof(T)*ys.size(), &ys[0]);
            q.enqueueReadBuffer(buffers[2], CL_FALSE, 0, sizeof(T)*zs.size(), &zs[0]);
         }

         // copies buffers to dst buffers on device
         void copy_to_dev(const cl::CommandQueue &q,
                          Buffer3D &dst,
                          const size_t amount_to_copy) const noexcept
         {
            q.enqueueCopyBuffer(buffers[0], dst.x(), 0, 0, amount_to_copy);
            q.enqueueCopyBuffer(buffers[1], dst.y(), 0, 0, amount_to_copy);
            q.enqueueCopyBuffer(buffers[2], dst.z(), 0, 0, amount_to_copy);
         }

         void zero_buffers(void)
         {
#ifdef CL_API_SUFFIX__VERSION_1_2
            const vcl_real_t zero = 0.0;
            for (const &buffer : buffers)
               vcl::queue.enqueueFillBuffer(buffer, &zero, sizeof(zero), buffer_size);
#else
            static cl::Kernel  zero_kernel;
            static bool initialised = false;

            if (!initialised)
            {
               zero_kernel = vcl::build_kernel_from_file("src/opencl/cl/zero_kernel.cl",
                                                      "zero_buffers",
                                                      vcl::context,
                                                      vcl::default_device);
               vcl::set_kernel_args(zero_kernel,
                                    buffers[0],
                                    buffers[1],
                                    buffers[2]);
               initialised = true;
            }

            vcl::kernel_call(zero_kernel, vcl::queue, buffer_size, vcl::local);

            vcl::queue.finish();
#endif // CL_API_SUFFIX__VERSION_1_2
         }

         // access each buffer, e.g. for passing to kernel
         cl::Buffer &x(void) noexcept { return buffers[0]; }
         cl::Buffer &y(void) noexcept { return buffers[1]; }
         cl::Buffer &z(void) noexcept { return buffers[2]; }

         // release memory by swapping with empty vector
         void free(void) noexcept
         {
            std::vector<cl::Buffer>().swap(buffers);
         }
      };

   }
}

#endif // VOPENCL_TYPEDEFS_HPP_
