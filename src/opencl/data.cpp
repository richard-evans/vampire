//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

#ifdef OPENCL_DEBUG
#include <fstream>
#endif

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
      cl::Context context;
      cl::CommandQueue queue;

#ifdef OPENCL_DEBUG
      std::ofstream OCLLOG("OpenCL.log");
#endif

      namespace rng
      {
         cl::Buffer urands;
         cl::Buffer grands;
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

      namespace stats
      {
         bool use_cpu = false;
         long counter(0L);

         cl::Buffer system_mask;
         cl::Buffer system_magnetization;
         cl::Buffer system_mean_magnetization;
         int system_mask_size(0);

         cl::Buffer material_mask;
         cl::Buffer material_magnetization;
         cl::Buffer material_mean_magnetization;
         int material_mask_size(0);

         cl::Buffer height_mask;
         cl::Buffer height_magnetization;
         cl::Buffer height_mean_magnetization;
         int height_mask_size(0);

         cl::Buffer material_height_mask;
         cl::Buffer material_height_magnetization;
         cl::Buffer material_height_mean_magnetization;
         int material_height_mask_size(0);
      }
   } // end of internal namespace
} // end of vopencl namespace

#endif // OPENCL
