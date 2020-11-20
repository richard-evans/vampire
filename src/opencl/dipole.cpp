//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <sstream>

#include "atoms.hpp"
#include "cells.hpp"
#include "demag.hpp"
#include "sim.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_utils.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      cl::Kernel update_dip;
      cl::Kernel update_atm_dip;
      cl::Kernel update_cell_mag;

      void update_dipolar_fields(void) noexcept
      {
         // dipole calculations enabled?
         if (!::dipole::activated) return;

         // check for previous demag update at same time
         if (::sim::time == ::demag::update_time) return;

         if (::sim::time % ::demag::update_rate != 0) return;

         ::demag::update_time = ::sim::time;

         // update cell magnetizations
         vcl::cells::mag_array.zero_buffer();
         vcl::kernel_call(update_cell_mag);
         vcl::queue.finish();

         // update cell dipolar fields
         vcl::kernel_call(update_dip);
         vcl::queue.finish();

         // update atomistic dipolar fields
         vcl::kernel_call(update_atm_dip);
      }
   }
}

#endif // OPENCL
