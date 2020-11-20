//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_TYPEDEFS_HPP_
#define VOPENCL_TYPEDEFS_HPP_

#include <cassert>
#include <type_traits>
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
      typedef cl_double  real_t;
      typedef cl_double3 real_t3;
#else
      typedef cl_float  real_t;
      typedef cl_float3 real_t3;
#endif // OPENCL_DP


      // Class for device buffer where each element is a 3 component vector.
      // i.e. stored as x0,y0,z0,x1,y1,z1,...
      // As the hardware will read from memory in chunks (typically 128bits)
      // this will allow a work item to read the x,y,z components in fewer
      // reads than if they were in different locations.

      class Buffer3D
      {
         std::vector<cl::Buffer> buff_container;

         // number of elements in one dimension
         size_t n_elems;

         size_t buffer_size;

#ifdef OPENCL_USE_VECTOR_TYPE
         typedef vcl::real_t3 elem_t;
         const static unsigned n_elems_per_set = 1;
#else
         typedef vcl::real_t elem_t;
         const static unsigned n_elems_per_set = 3;
#endif

      public:

         Buffer3D(void) {}

         // initialize without writing, but with size
         // e.g. for use when generating buffer on device
         Buffer3D(const cl::Context &c,
                  const cl_mem_flags fs,
                  const size_t n) noexcept
            : buff_container(1), n_elems(n), buffer_size(n * n_elems_per_set * sizeof (elem_t))
         {
            buff_container[0] = cl::Buffer(c, fs, buffer_size);
         }

         // initialize with 3 vectors to write data to device
         // use template so that class can be initialised
         // with std::vector<double> or std::vector<float>
         template <typename R>
         Buffer3D(const cl::Context &c,
                  const cl::CommandQueue &q,
                  const cl_mem_flags fs,
                  const std::vector<R> &xs,
                  const std::vector<R> &ys,
                  const std::vector<R> &zs) noexcept
            : buff_container(1)
         {
            assert(xs.data() != nullptr   &&
                   xs.size() == ys.size() &&
                   xs.size() == zs.size());

            n_elems = xs.size();
            buffer_size = n_elems * n_elems_per_set * sizeof(elem_t);

            std::vector<elem_t> buff(n_elems_per_set * n_elems);
            for (size_t i=0; i<n_elems; ++i)
            {
#ifdef OPENCL_USE_VECTOR_TYPE
               buff[i] = elem_t{xs[i], ys[i], zs[i]};
#else
               buff[3*i+0] = vcl::real_t(xs[i]);
               buff[3*i+1] = vcl::real_t(ys[i]);
               buff[3*i+2] = vcl::real_t(zs[i]);
#endif // OPENCL_USE_VECTOR_TYPE
            }

            buff_container[0] = vcl::create_device_buffer(buff, fs, CL_TRUE, q, c);
         }

         // reads data from device, assumes host vectors already have enough capacity
         template <typename R>
         void copy_to_host(const cl::CommandQueue &q,
                           std::vector<R> &xs,
                           std::vector<R> &ys,
                           std::vector<R> &zs) const noexcept
         {
            if (n_elems == 0) return;

            assert(xs.size() == n_elems &&
                   ys.size() == n_elems &&
                   zs.size() == n_elems);

            std::vector<elem_t> buff(n_elems_per_set * n_elems);
            q.enqueueReadBuffer(buff_container[0], CL_TRUE, 0, buffer_size, buff.data());

            for (size_t i=0; i<n_elems; ++i)
            {
#ifdef OPENCL_USE_VECTOR_TYPE
               xs[i] = R(buff[i].s[0]);
               ys[i] = R(buff[i].s[1]);
               zs[i] = R(buff[i].s[2]);
#else
               xs[i] = R(buff[3*i+0]);
               ys[i] = R(buff[3*i+1]);
               zs[i] = R(buff[3*i+2]);
#endif // OPENCL_USE_VECTOR_TYPE
            }
         }

         // copies buffer to dst buffer on device
         void copy_to_dev(const cl::CommandQueue &q,
                          Buffer3D &dst) const noexcept
         {
            q.enqueueCopyBuffer(buff_container[0], dst.buff_container[0], 0, 0, buffer_size);
            q.finish();
         }

         // writes zeros over the data
         void zero_buffer(void) noexcept
         {
#ifdef CL_API_SUFFIX__VERSION_1_2
            const elem_t zero{0.0};
            vcl::queue.enqueueFillBuffer(buff_container[0], &zero, sizeof (elem_t), n_elems_per_set * n_elems);
#else
            const std::vector<elem_t> zeros(3 * n_elems, {0.0});
            vcl::queue.enqueueWriteBuffer(buff_container[0], CL_FALSE, 0, buffer_size, zeros.data());
#endif // CL_API_SUFFIX__VERSION_1_2

            vcl::queue.finish();
         }

         // access the buffer (to pass to kernels)
         cl::Buffer &buffer(void) noexcept { return buff_container[0]; }

         // release memory by swapping with empty vector
         void free(void) noexcept
         {
            std::vector<cl::Buffer>().swap(buff_container);
         }
      };

   }
}

#endif // VOPENCL_TYPEDEFS_HPP_
