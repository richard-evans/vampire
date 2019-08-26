//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

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
