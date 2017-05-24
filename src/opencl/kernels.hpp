//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef KERNELS_HPP_
#define KERNELS_HPP_

#include "opencl_include.hpp"

namespace vopencl
{
   namespace internal
   {
      namespace llg
      {
         extern cl::Kernel predictor_step;
         extern cl::Kernel corrector_step;
      }

      // dipole
      extern cl::Kernel update_dip;
      extern cl::Kernel update_atm_dip;
      extern cl::Kernel update_cell_mag;

      // external fields
      extern cl::Kernel update_ext;

      // spin fields
      extern cl::Kernel update_nexch_spin_fields;

      // exchange
      namespace exchange
      {
         extern cl::Kernel calculate_exchange;
      }

      // gaussian rng
      namespace rng
      {
         extern cl::Kernel grng;
      }
   }
}

#endif // KERNELS_HPP_
