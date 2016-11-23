//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------
//

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
         extern cl::Buffer urands;
         extern cl::Buffer grands;
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

      namespace stats
      {
         extern bool use_cpu;
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

   } // end of internal namespace
#endif // OPENCL

} // end of vopencl namespace

#endif // VOPENCL_DATA_HPP_
