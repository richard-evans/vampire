#ifndef VOPENCL_DEVICE_VECTOR_HPP_
#define VOPENCL_DEVICE_VECTOR_HPP_

#include <vector>
#include "opencl_include.hpp"

namespace vopencl
{
   template <typename T>
   class DeviceVector
   {
      cl::Context context;
      cl_mem_flags flags;

   public:
      cl::CommandQueue queue;
      cl::Buffer buffer;
      size_t buffer_size;

      DeviceVector(const cl::Context &_context,
                   const cl::CommandQueue &_queue,
                   const cl_mem_flags _flags,
                   const size_t N) :
         context(_context),
         flags(_flags),
         queue(_queue),
         buffer_size(N * sizeof(T))
      {
         cl::Buffer buffer = cl::Buffer(context, flags, buffer_size);
      }

      void operator=(std::vector<T> &HV)
      {
         queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, buffer_size, &HV[0]);
      }
   };
}

#endif // VOPENCL_DEVICE_VECTOR_HPP_
