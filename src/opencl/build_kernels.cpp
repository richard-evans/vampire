//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <sstream>
#include <string>

#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "sim.hpp"

#include "internal.hpp"
#include "kernels.hpp"
#include "opencl_utils.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      void build_kernels(void) noexcept
      {
         std::string default_opts("-Isrc/opencl/cl "
                                  "-cl-fast-relaxed-math "
#ifdef OPENCL_DP
                                  "-DOPENCL_DP "
#endif
#ifdef OPENCL_USE_NATIVE_FUNCTIONS
                                  "-DOPENCL_USE_NATIVE_FUNCTIONS "
#endif
#ifdef OPENCL_USE_VECTOR_TYPE
                                  "-DOPENCL_USE_VECTOR_TYPE "
#endif
            );
         std::ostringstream params;
         params << "-DNUM_ATOMS=" << ::atoms::num_atoms << ' ';
         params << "-DNUM_CELLS=" << ::cells::num_cells << ' ';
         params << "-DDT=" << ::mp::dt << ' ';

         default_opts.append(params.str());

         // dipole calculations enabled?
         if (::dipole::activated)
         {
            // build dipolar kernels
            vcl::update_dip =      vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                               "update_dipole_fields",
                                                               default_opts);
            vcl::update_atm_dip =  vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                               "update_atm_dipole_fields",
                                                               default_opts);
            vcl::update_cell_mag = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                               "update_cell_magnetization",
                                                               default_opts);
         }

         // external field kernel
         vcl::update_ext = vcl::build_kernel_from_file("src/opencl/cl/external_fields.cl",
                                                       "update_external_fields",
                                                       default_opts);

         // matrix multiplication for exchange
         std::ostringstream exch_opts;
         exch_opts << default_opts << "-DEXCH_TYPE=" << ::atoms::exchange_type;
         vcl::exchange::calculate_exchange = vcl::build_kernel_from_file("src/opencl/cl/csrmatmul.cl",
                                                                         "matmul",
                                                                         exch_opts.str());

         // non exchange spin fields
         vcl::update_nexch_spin_fields = vcl::build_kernel_from_file("src/opencl/cl/spin_fields.cl",
                                                                     "update_nexch_spin_fields",
                                                                     default_opts);

         // gaussian prng
         vcl::rng::grng = vcl::build_kernel_from_file("src/opencl/cl/random.cl",
                                                      "gen_grands",
                                                      default_opts
#ifdef ENABLE_MULTIPLE_DEVICES
                                                      ,
                                                      vcl::extra_device,
                                                      vcl::context_other
#endif // ENABLE_MULTIPLE_DEVICES
            );

         // llg heun steps
         vcl::llg::predictor_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                                "llg_heun_predictor_step",
                                                                default_opts);
         vcl::llg::corrector_step = vcl::build_kernel_from_file("src/opencl/cl/llg_heun.cl",
                                                                "llg_heun_corrector_step",
                                                                default_opts);
      }
   }
}
#endif // OPENCL
