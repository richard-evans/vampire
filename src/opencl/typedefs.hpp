#ifndef VOPENCL_TYPEDEFS_HPP_
#define VOPENCL_TYPEDEFS_HPP_

#include <vector>

#include "opencl_include.hpp"

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
         std::vector<cl::Buffer> buf;

      public:

         Buffer3D(void) : buf(3) {}

         // initialize without writing, but with size
         // e.g. for use when generating buffers on device
         Buffer3D(cl::Context &c,
                  cl_mem_flags fs,
                  size_t size)
            : buf(3)
         {
            buf[0] = cl::Buffer(c, fs, size);
            buf[1] = cl::Buffer(c, fs, size);
            buf[2] = cl::Buffer(c, fs, size);
         }
            
         // initialize with 3 vectors to write data to device
         template <typename T>
         Buffer3D(cl::Context &c,
                  cl::CommandQueue &q,
                  cl_mem_flags fs,
                  std::vector<T> &xs,
                  std::vector<T> &ys,
                  std::vector<T> &zs)
            : buf(3)
         {
            buf[0] = cl::Buffer(c, fs, sizeof(T)*xs.size());
            buf[1] = cl::Buffer(c, fs, sizeof(T)*ys.size());
            buf[2] = cl::Buffer(c, fs, sizeof(T)*zs.size());

            q.enqueueWriteBuffer(buf[0], CL_FALSE, 0, sizeof(T)*xs.size(), &xs[0]);
            q.enqueueWriteBuffer(buf[1], CL_FALSE, 0, sizeof(T)*ys.size(), &ys[0]);
            q.enqueueWriteBuffer(buf[2], CL_FALSE, 0, sizeof(T)*zs.size(), &zs[0]);
         }

         // reads data from device, assumes host vectors already have enough capacity
         template <typename T>
         void copy_to_host(cl::CommandQueue &q,
                               std::vector<T> &xs,
                               std::vector<T> &ys,
                               std::vector<T> &zs)
         {
            q.enqueueReadBuffer(buf[0], CL_FALSE, 0, sizeof(T)*xs.size(), &xs[0]);
            q.enqueueReadBuffer(buf[1], CL_FALSE, 0, sizeof(T)*ys.size(), &ys[0]);
            q.enqueueReadBuffer(buf[2], CL_FALSE, 0, sizeof(T)*zs.size(), &zs[0]);
         }

         // copies buffers to dst buffers on device
         void copy_to_dev(cl::CommandQueue &q,
                          Buffer3D &dst,
                          size_t amount_to_copy)
         {
            q.enqueueCopyBuffer(buf[0], dst.x(), 0, 0, amount_to_copy);
            q.enqueueCopyBuffer(buf[1], dst.y(), 0, 0, amount_to_copy);
            q.enqueueCopyBuffer(buf[2], dst.z(), 0, 0, amount_to_copy);
         }

         // access each buffer, e.g. for passing to kernel
         cl::Buffer &x(void) { return buf[0]; }
         cl::Buffer &y(void) { return buf[1]; }
         cl::Buffer &z(void) { return buf[2]; }

         // release memory by swapping with empty vector
         void free(void)
         {
            std::vector<cl::Buffer>().swap(buf);
         }
      };

   }
}

#endif // VOPENCL_TYPEDEFS_HPP_
