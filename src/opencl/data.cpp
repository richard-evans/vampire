//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

#ifdef OPENCL_DEBUG
#include <fstream>
#endif // OPENCL_DEBUG

// Vampire headers
#include "vopencl.hpp"

// vopencl module headers
#include "data.hpp"
#include "internal.hpp"
#include "opencl_include.hpp"
#include "statistics.hpp"
#include "typedefs.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   //------------------------------------------------------------------------
   // Shared variables inside vopencl module
   //------------------------------------------------------------------------

   namespace internal
   {
      cl::Device default_device;
      cl::Device extra_device;
      cl::Context context;
      cl::Context context_other;
      cl::CommandQueue queue;
      cl::CommandQueue queue_other;

#ifdef OPENCL_LOG
      std::ofstream OCLLOG("OpenCL.log");
#endif // OPENCL_LOG

      namespace rng
      {
         cl::Buffer state;
         cl::Buffer grands;
         cl::Buffer grands_copy;

         unsigned n_rands;
      }

      namespace atoms
      {
         vcl::Buffer3D spin_array;

         vcl::Buffer3D coord_array;

         cl::Buffer type_array;

         cl::Buffer cell_array;

         cl::Buffer limits;
         cl::Buffer neighbours;

         cl::Buffer spin_norm_array;
      }

      namespace cells
      {
         vcl::Buffer3D coord_array;

         vcl::Buffer3D mag_array;

         vcl::Buffer3D field_array;

         cl::Buffer volume_array;

         cl::Buffer num_atoms;
      }

      namespace mp
      {
         cl::Buffer materials;
      }

      vcl::Buffer3D total_spin_field_array;

      vcl::Buffer3D total_external_field_array;

      vcl::Buffer3D dipolar_field_array;

      namespace llg
      {
         vcl::Buffer3D spin_buffer_array;
         vcl::Buffer3D dS_array;
         cl::Buffer heun_parameters_device;
      }

      namespace exchange
      {
         cl::Buffer J_vals_d;
         cl::Buffer Jxx_vals_d;
         cl::Buffer Jyy_vals_d;
         cl::Buffer Jzz_vals_d;
      }

      namespace stats
      {
         bool use_cpu = true;
      }
   } // end of internal namespace
} // end of vopencl namespace

#endif // OPENCL
