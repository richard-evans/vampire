#include <sstream>

#include "atoms.hpp"
#include "cells.hpp"
#include "demag.hpp"
#include "sim.hpp"

#include "data.hpp"
#include "internal.hpp"
#include "opencl_utils.hpp"

#ifdef OPENCL

namespace vcl = ::vopencl::internal;

namespace vopencl
{
   namespace internal
   {
      bool compiled_update_dip = false;
      bool compiled_update_atm_dip = false;
      cl::Kernel update_dip;
      cl::Kernel update_atm_dip;

      void update_dipolar_fields()
      {
         // dipole calculations enabled?
         if (::sim::hamiltonian_simulation_flags[4]!=1) return;

         // check for previous demag update at same time
         if (::sim::time == ::demag::update_time) return;

         if (::sim::time % ::demag::update_rate != 0) return;

         ::demag::update_time = ::sim::time;

         //update_cell_magnetizations();

         if (!compiled_update_dip)
         {
            std::ostringstream opts;
            opts << "-DN=" << ::cells::num_cells;
            update_dip = vcl::build_kernel_from_file("dipole.cl",
                                                     "update_dipole_fields",
                                                     vcl::context, vcl::default_device,
                                                     opts.str());
            compiled_update_dip = true;
         }

         if (!compiled_update_atm_dip)
         {
            std::ostringstream opts;
            opts << "-DN=" << ::atoms::num_atoms;
            update_atm_dip = vcl::build_kernel_from_file("dipole.cl",
                                                         "update_atm_dipole_fields",
                                                         vcl::context, vcl::default_device,
                                                         opts.str());
            compiled_update_atm_dip = true;
         }

         cl::CommandQueue update_q(vcl::context, vcl::default_device);

         cl::NDRange global(::atoms::num_atoms);
         cl::NDRange local(0);

         // update cell dipolar fields
         vcl::kernel_call(update_dip, update_q, global, local,
                          vcl::cells::x_mag_array,
                          vcl::cells::y_mag_array,
                          vcl::cells::z_mag_array,
                          vcl::cells::x_coord_array,
                          vcl::cells::y_coord_array,
                          vcl::cells::z_coord_array,
                          vcl::cells::volume_array,
                          vcl::cells::x_field_array,
                          vcl::cells::y_field_array,
                          vcl::cells::z_field_array);

         // update atomistic dipolar fields
         vcl::kernel_call(update_atm_dip, update_q, global, local,
                          vcl::cells::x_field_array,
                          vcl::cells::y_field_array,
                          vcl::cells::z_field_array,
                          vcl::x_dipolar_field_array,
                          vcl::y_dipolar_field_array,
                          vcl::z_dipolar_field_array,
                          vcl::atoms::cell_array);

         update_q.finish();
      }
   }
}

#endif // OPENCL
