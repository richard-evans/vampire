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
#include "opencl_include.hpp"
#include "vopencl.hpp"

// vopencl module headers
#include "internal.hpp"
#include "data.hpp"
#include "statistics.hpp"

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
      cl::Device default_device;
      cl::Context context;

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
         cl::Buffer x_spin_array;
         cl::Buffer y_spin_array;
         cl::Buffer z_spin_array;

         cl::Buffer x_coord_array;
         cl::Buffer y_coord_array;
         cl::Buffer z_coord_array;

         cl::Buffer type_array;

         cl::Buffer cell_array;

         cl::Buffer limits;
         cl::Buffer neighbours;

         cl::Buffer spin_norm_array;
      }

      namespace cells
      {
         cl::Buffer x_coord_array;
         cl::Buffer y_coord_array;
         cl::Buffer z_coord_array;

         cl::Buffer x_mag_array;
         cl::Buffer y_mag_array;
         cl::Buffer z_mag_array;

         cl::Buffer x_field_array;
         cl::Buffer y_field_array;
         cl::Buffer z_field_array;

         cl::Buffer volume_array;

         cl::Buffer num_atoms;
      }

      namespace mp
      {
         cl::Buffer materials;
      }

      cl::Buffer x_total_spin_field_array;
      cl::Buffer y_total_spin_field_array;
      cl::Buffer z_total_spin_field_array;

      cl::Buffer x_total_external_field_array;
      cl::Buffer y_total_external_field_array;
      cl::Buffer z_total_external_field_array;

      cl::Buffer x_dipolar_field_array;
      cl::Buffer y_dipolar_field_array;
      cl::Buffer z_dipolar_field_array;

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
#endif // OPENCL

} // end of vopencl namespace

