//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "atoms.hpp"
#include "vopencl.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace config
   {
      void synchronise(void)
      {
         vcl::atoms::spin_array.copy_to_host(vcl::queue,
                                             ::atoms::x_spin_array,
                                             ::atoms::y_spin_array,
                                             ::atoms::z_spin_array);
      }
   }
}

#endif
