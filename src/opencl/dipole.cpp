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

      void update_dipolar_fields()
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
            update_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                     "update_dipole_fields",
                                                     vcl::context, vcl::default_device,
                                                     opts.str());
            compiled_update_dip = true;
         }

         if (!compiled_update_atm_dip)
         {
            std::ostringstream opts;
            opts << "-DN_ATOMS=" << ::atoms::num_atoms;
            update_atm_dip = vcl::build_kernel_from_file("src/opencl/cl/dipole.cl",
                                                         "update_atm_dipole_fields",
                                                         vcl::context, vcl::default_device,
                                                         opts.str());
            compiled_update_atm_dip = true;
         }

         cl::CommandQueue update_q(vcl::context, vcl::default_device);

         cl::NDRange global(::atoms::num_atoms);

         // update cell dipolar fields
         vcl::kernel_call(update_dip, update_q, global, vcl::local,
                          vcl::cells::mag_array.x(),
                          vcl::cells::mag_array.y(),
                          vcl::cells::mag_array.z(),
                          vcl::cells::coord_array.x(),
                          vcl::cells::coord_array.y(),
                          vcl::cells::coord_array.z(),
                          vcl::cells::volume_array,
                          vcl::cells::field_array.x(),
                          vcl::cells::field_array.y(),
                          vcl::cells::field_array.z());

         update_q.finish();

         // update atomistic dipolar fields
         vcl::kernel_call(update_atm_dip, update_q, global, vcl::local,
                          vcl::cells::field_array.x(),
                          vcl::cells::field_array.y(),
                          vcl::cells::field_array.z(),
                          vcl::dipolar_field_array.x(),
                          vcl::dipolar_field_array.y(),
                          vcl::dipolar_field_array.z(),
                          vcl::atoms::cell_array);

         update_q.finish();
      }

      void update_cell_magnetizations(void)
      {
         cl::CommandQueue cell_q(vcl::context, vcl::default_device);
         cl::NDRange global(::cells::num_cells);

         size_t buff_size = ::cells::num_cells * sizeof(vcl_real_t);
         vcl_real_t zero = 0;
         cell_q.enqueueFillBuffer(vcl::cells::mag_array.x(), &zero, sizeof(zero), buff_size);
         cell_q.enqueueFillBuffer(vcl::cells::mag_array.y(), &zero, sizeof(zero), buff_size);
         cell_q.enqueueFillBuffer(vcl::cells::mag_array.z(), &zero, sizeof(zero), buff_size);

         if (!compiled_update_cell_magnetization)
         {
            std::ostringstream opts;
            opts << "-DN_ATOMS=" << ::atoms::num_atoms;
            update_cell_mag = vcl::build_kernel_from_file("dipole.cl",
                                                          "update_cell_magnetization",
                                                          vcl::context, vcl::default_device,
                                                          opts.str());
            compiled_update_cell_magnetization = true;
         }
         
         cell_q.finish();

         vcl::kernel_call(update_cell_mag, cell_q, global, vcl::local,
                          vcl::atoms::spin_array.x(),
                          vcl::atoms::spin_array.y(),
                          vcl::atoms::spin_array.z(),
                          vcl::atoms::type_array,
                          vcl::atoms::cell_array,
                          vcl::mp::materials,
                          vcl::cells::mag_array.x(),
                          vcl::cells::mag_array.y(),
                          vcl::cells::mag_array.z());
      }
   }
}

#endif // OPENCL
