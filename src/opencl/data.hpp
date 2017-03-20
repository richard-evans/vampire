//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#ifndef VOPENCL_DATA_HPP_
#define VOPENCL_DATA_HPP_

// C++ standard library headers

// Vampire headers
#include "vopencl.hpp"

// vopencl module headers
#include "internal.hpp"
#include "opencl_include.hpp"
#include "typedefs.hpp"

namespace vopencl
{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   //------------------------------------------------------------------------
   // Shared variables inside vopencl module
   //------------------------------------------------------------------------


#ifdef OPENCL
   namespace internal
   {
      namespace rng
      {
         extern cl::Buffer state;
         extern cl::Buffer grands;
         extern cl::Buffer grands_copy;

         extern unsigned n_rands;
      }

      namespace atoms
      {
         extern Buffer3D spin_array;

         extern Buffer3D coord_array;

         extern cl::Buffer type_array;

         extern cl::Buffer cell_array;

         extern cl::Buffer limits;
         extern cl::Buffer neighbours;

         extern cl::Buffer spin_norm_array;
      }

      namespace cells
      {
         extern Buffer3D coord_array;

         extern Buffer3D mag_array;

         extern Buffer3D field_array;

         extern cl::Buffer volume_array;

         extern cl::Buffer num_atoms;
      }

      namespace mp
      {
         extern cl::Buffer materials;
      }

      extern Buffer3D total_spin_field_array;

      extern Buffer3D total_external_field_array;

      extern Buffer3D dipolar_field_array;

      namespace llg
      {
         extern Buffer3D spin_buffer_array;
         extern Buffer3D dS_array;
         extern cl::Buffer heun_parameters_device;
      }

      namespace exchange
      {
         extern cl::Buffer J_vals_d;
         extern cl::Buffer Jxx_vals_d;
         extern cl::Buffer Jyy_vals_d;
         extern cl::Buffer Jzz_vals_d;
      }

      namespace stats
      {
         extern bool use_cpu;
      }

   } // end of internal namespace
#endif // OPENCL

} // end of vopencl namespace

#endif // VOPENCL_DATA_HPP_
