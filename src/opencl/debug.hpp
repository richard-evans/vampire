#ifndef VOPENCL_DEBUG_HPP_
#define VOPENCL_DEBUG_HPP_

#include <cassert>
#include <cstring>
#include <vector>

#include "internal.hpp"
#include "opencl_include.hpp"

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      namespace debug
      {
         template <typename T>
         void verify_copy(const cl::Buffer     &dev_vector,
                          const std::vector<T> &host_vector)
         {
            const size_t n_elems = host_vector.size();
            if (n_elems == 0)
               return;

            const size_t buffer_size = n_elems * sizeof (T);
            std::vector<T> copy_vector(n_elems);

            vcl::queue.enqueueReadBuffer(dev_vector, CL_TRUE, 0, buffer_size, copy_vector.data());

            assert(std::memcmp(host_vector.data(), copy_vector.data(), buffer_size) == 0);
         }
      }
   }
}

#endif // VOPENCL_DEBUG_HPP_
