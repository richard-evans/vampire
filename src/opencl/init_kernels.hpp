#ifndef INIT_KERNELS_HPP_
#define INIT_KERNELS_HPP_

namespace vopencl
{
   namespace internal
   {
      void build_kernels(void) noexcept;
      bool initialize_kernels(void) noexcept;
   }
}

#endif // INIT_KERNELS_HPP_
