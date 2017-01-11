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
      bool compiled_update_dip = false;
      bool compiled_update_atm_dip = false;
      bool compiled_update_cell_magnetization = false;
      cl::Kernel update_dip;
      cl::Kernel update_atm_dip;
      cl::Kernel update_cell_mag;

      void update_dipolar_fields(void) noexcept
      {
         // dipole calculations enabled?
         if (::sim::hamiltonian_simulation_flags[4]!=1) return;

         // check for previous demag update at same time
         if (::sim::time == ::demag::update_time) return;

         if (::sim::time % ::demag::update_rate != 0) return;

         ::demag::update_time = ::sim::time;

         update_cell_magnetizations();

         if (!compiled_update_dip)
         {
            std::ostringstream opts;
            opts << "-DN_CELLS=" << ::cells::num_cells;
            opts << " -DN_ATOMS=" << ::atoms::num_atoms;
            update_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                     "update_dipole_fields",
                                                     vcl::context, vcl::default_device,
                                                     opts.str());
            compiled_update_dip = true;

            vcl::set_kernel_args(update_dip,
                                 vcl::cells::mag_array.buffer(),
                                 vcl::cells::coord_array.buffer(),
                                 vcl::cells::volume_array,
                                 vcl::cells::field_array.buffer());
         }

         if (!compiled_update_atm_dip)
         {
            std::ostringstream opts;
            opts << "-DN_CELLS=" << ::cells::num_cells;
            opts << " -DN_ATOMS=" << ::atoms::num_atoms;
            update_atm_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                         "update_atm_dipole_fields",
                                                         vcl::context, vcl::default_device,
                                                         opts.str());
            compiled_update_atm_dip = true;

            vcl::set_kernel_args(update_atm_dip,
                                 vcl::cells::field_array.buffer(),
                                 vcl::dipolar_field_array.buffer(),
                                 vcl::atoms::cell_array);
         }

         const cl::NDRange global(::atoms::num_atoms);

         // update cell dipolar fields
         vcl::kernel_call(update_dip, vcl::queue, global, vcl::local);

         vcl::queue.finish();

         // update atomistic dipolar fields
         vcl::kernel_call(update_atm_dip, vcl::queue, global, vcl::local);

         vcl::queue.finish();

      }

      void update_cell_magnetizations(void) noexcept
      {
         const cl::CommandQueue cell_q(vcl::context, vcl::default_device);
         const cl::NDRange global(::cells::num_cells);

         vcl::cells::mag_array.zero_buffer();

         if (!compiled_update_cell_magnetization)
         {
            std::ostringstream opts;
            opts << "-DN_ATOMS=" << ::atoms::num_atoms;
            update_cell_mag = vcl::build_kernel_from_file("dipole.cl",
                                                          "update_cell_magnetization",
                                                          vcl::context, vcl::default_device,
                                                          opts.str());
            compiled_update_cell_magnetization = true;

            vcl::set_kernel_args(update_cell_mag,
                                 vcl::atoms::spin_array.buffer(),
                                 vcl::atoms::type_array,
                                 vcl::atoms::cell_array,
                                 vcl::mp::materials,
                                 vcl::cells::mag_array.buffer());
         }

         vcl::kernel_call(update_cell_mag, cell_q, global, vcl::local);

      }
   }
}

#endif // OPENCL
