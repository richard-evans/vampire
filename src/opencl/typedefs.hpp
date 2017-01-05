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
      typedef cl_double  vcl_real_t;
      typedef cl_double3 vcl_real_vec_t;
#else
      typedef cl_float  vcl_real_t;
      typedef cl_float3 vcl_real_vec_t;
#endif // OPENCL_DP


      // Class for set of three buffers, i.e. for x,y,z arrays
      // uses non blocking reads/writes where possible

      // T is standard C++ type (i.e. float)
      // T_vec is OpenCL type (i.e. float3)
      template <typename T, typename T_vec>
      class Buffer3D
      {
         std::vector<cl::Buffer> buffer;
         size_t n_elems;

      public:

         Buffer3D(void) noexcept : buffer(1), n_elems(0) {}

         // initialize without writing, but with size
         // e.g. for use when generating buffer on device
         Buffer3D(const cl::Context &c,
                  const cl_mem_flags fs,
                  const size_t n) noexcept
            : buffer(1), n_elems(n)
         {
            buffer[0] = cl::Buffer(c, fs, n*sizeof(T_vec));
         }

          // initialize with 3 vectors to write data to device
         template <typename R>
         Buffer3D(const cl::Context &c,
                  const cl::CommandQueue &q,
                  const cl_mem_flags fs,
                  const std::vector<R> &xs,
                  const std::vector<R> &ys,
                  const std::vector<R> &zs) noexcept
            : buffer(1)
         {
            n_elems = xs.size();

            buffer[0] = cl::Buffer(c, fs, sizeof(T_vec)*xs.size());

            std::vector<T_vec> buff(n_elems);
            for (size_t i=0; i<n_elems; ++i)
            {
               buff[i] = {T(xs[i]), T(ys[i]), T(zs[i])};
            }

            q.enqueueWriteBuffer(buffer[0], CL_FALSE, 0, sizeof(T_vec)*n_elems, &buff[0]);
         }

         // reads data from device, assumes host vectors already have enough capacity
         template <typename R>
         void copy_to_host(const cl::CommandQueue &q,
                           std::vector<R> &xs,
                           std::vector<R> &ys,
                           std::vector<R> &zs) const noexcept
         {
            std::vector<T_vec> buff(n_elems);
            q.enqueueReadBuffer(buffer[0], CL_FALSE, 0, sizeof(T_vec)*n_elems, &buff[0]);

            q.finish();

            for (size_t i=0; i<n_elems; ++i)
            {
               xs[i] = R(buff[i].s[0]);
               ys[i] = R(buff[i].s[1]);
               zs[i] = R(buff[i].s[2]);
            }
         }

         // copies buffer to dst buffer on device
         void copy_to_dev(const cl::CommandQueue &q,
                          Buffer3D &dst,
                          const size_t amount_to_copy) const noexcept
         {
            q.enqueueCopyBuffer(buffer[0], dst.buffer[0], 0, 0, amount_to_copy);
         }

         void zero_buffer(void) noexcept
         {
#ifdef CL_API_SUFFIX__VERSION_1_2
            const T_vec zero = {0};
            vcl::queue.enqueueFillBuffer(buffer[0], &zero, sizeof T_vec, n_elems);
#else
            const std::vector<T_vec> zeros(n_elems, {0});
            vcl::queue.enqueueWriteBuffer(buffer[0], CL_FALSE, 0, sizeof(T_vec) * n_elems, &zeros[0]);
#endif // CL_API_SUFFIX__VERSION_1_2

            vcl::queue.finish();
         }

         // release memory by swapping with empty vector
         void free(void) noexcept
         {
            std::vector<cl::Buffer>().swap(buffer);
         }
      };

   }
}

#endif // VOPENCL_TYPEDEFS_HPP_
