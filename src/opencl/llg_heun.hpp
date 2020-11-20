//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_LLG_HEUN_HPP_
#define VOPENCL_LLG_HEUN_HPP_

#include "internal.hpp"
#include "opencl_include.hpp"
#include "typedefs.hpp"

namespace vopencl
{
   namespace internal
   {
      struct heun_parameter_t
      {
         vcl::real_t prefactor;
         vcl::real_t lambda_times_prefactor;
      };

      namespace llg
      {
         extern bool initialized;

         extern cl::Buffer x_spin_array;
         extern cl::Buffer y_spin_array;
         extern cl::Buffer z_spin_array;

         extern cl::Buffer dS_x_array;
         extern cl::Buffer dS_y_array;
         extern cl::Buffer dS_z_array;

         extern cl::Buffer heun_parameters_device;

         void init(void) noexcept;
         void step(void) noexcept;
      }
   }
}

#endif // VOPENCL_LLG_HEUN_HPP_
