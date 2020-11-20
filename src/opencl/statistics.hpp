//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_STATISTICS_HPP_
#define VOPENCL_STATISTICS_HPP_

#include "stats.hpp"

#include "opencl_include.hpp"
#include "internal.hpp"

namespace vopencl
{
#ifdef OPENCL
   namespace internal
   {
      namespace stats
      {
         extern long counter;

         extern cl::Buffer system_mask;
         extern cl::Buffer system_magnetization;
         extern cl::Buffer system_mean_magnetization;
         extern int system_mask_size;

         extern cl::Buffer material_mask;
         extern cl::Buffer material_magnetization;
         extern cl::Buffer material_mean_magnetization;
         extern int material_mask_size;

         extern cl::Buffer height_mask;
         extern cl::Buffer height_magnetization;
         extern cl::Buffer height_mean_magnetization;
         extern int height_mask_size;

         extern cl::Buffer material_height_mask;
         extern cl::Buffer material_height_magnetization;
         extern cl::Buffer material_height_mean_magnetization;
         extern int material_height_mask_size;
      }
   }
#endif // OPENCL
}

#endif // VOPENCL_STATISTICS_HPP_
