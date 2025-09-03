//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "../../hdr/cuda.hpp"

// Local cuda headers
#include "data.hpp"
#include "llg_heun.hpp"
#include "internal.hpp"

#include "exchange_fields.hpp"
#include "statistics.hpp"
#include "cuda_utils.hpp"

#include "monte_carlo.hpp"

#include <cuda_profiler_api.h>

#ifdef CUDA
namespace cu = vcuda::internal;
#endif

namespace vcuda{

   //----------------------------------------------------------------------------------
   // Function de-initialize gpu data
   //----------------------------------------------------------------------------------
   void finalize(){

      // Only compile code if CUDA enabled
#ifdef CUDA

      vcuda::internal::__finalize ();

      vcuda::internal::mc::finalise();
      cudaProfilerStop();
#endif

      return;
   }

#ifdef CUDA
   namespace internal
   {
      void __finalize ()
      {

         check_cuda_errors (__FILE__, __LINE__);
         std::cout << "CUDA time taken \t" << cuda_timer.seconds_elapsed() << std::endl;
         // De-allocate the exchange fields
         cu::exchange::finalise_exchange();

         check_cuda_errors (__FILE__, __LINE__);
         // Free all device memory

         // Free rand states
         cudaFree(d_rand_state);
         // Free atom topologies
         cudaFree(atoms::d_neighbours);
         cudaFree(atoms::d_limits);
         // Free material parameters
         cudaFree(mp::d_material_params_r);
         cudaFree(mp::d_material_params);
         // Free dipole variables
         cudaFree(cells::d_tensor_zz);
         cudaFree(cells::d_tensor_yz);
         cudaFree(cells::d_tensor_yy);
         cudaFree(cells::d_tensor_xz);
         cudaFree(cells::d_tensor_xy);
         cudaFree(cells::d_tensor_xx);

         cudaFree(cells::d_num_atoms_in_cell);

         cudaFree(cells::d_z_cell_mu0H_field);
         cudaFree(cells::d_y_cell_mu0H_field);
         cudaFree(cells::d_x_cell_mu0H_field);
         cudaFree(cells::d_z_cell_field);
         cudaFree(cells::d_y_cell_field);
         cudaFree(cells::d_x_cell_field);
         // Free cell variables
         cudaFree(cells::d_cell_id_array);
         cudaFree(cells::d_num_atoms);
         cudaFree(cells::d_volume);

         cudaFree(cells::d_z_mag);
         cudaFree(cells::d_y_mag);
         cudaFree(cells::d_x_mag);

         cudaFree(cells::d_z_coord);
         cudaFree(cells::d_y_coord);
         cudaFree(cells::d_x_coord);
         // Free neel tensor
         cudaFree(d_neel_tensor);
         // Free field variables
         cudaFree(d_z_mu0H_dip_field);
         cudaFree(d_y_mu0H_dip_field);
         cudaFree(d_x_mu0H_dip_field);
         cudaFree(d_z_dip_field);
         cudaFree(d_y_dip_field);
         cudaFree(d_x_dip_field);
         cudaFree(d_z_external_field);
         cudaFree(d_y_external_field);
         cudaFree(d_x_external_field);
         cudaFree(d_z_spin_field);
         cudaFree(d_y_spin_field);
         cudaFree(d_x_spin_field);
         cudaFree(d_spin_field);
         // Free atom variables
         cudaFree(atoms::d_cells);
         cudaFree(atoms::d_materials);
         // Free atom coords
         cudaFree(atoms::d_z_coord);
         cudaFree(atoms::d_y_coord);
         cudaFree(atoms::d_x_coord);
         // Free spin transfer buffers
         cudaFree(h_z_spin_transfer_buffer);
         cudaFree(h_y_spin_transfer_buffer);
         cudaFree(h_x_spin_transfer_buffer);
         // Free atom spins
         cudaFree(atoms::d_z_spin);
         cudaFree(atoms::d_y_spin);
         cudaFree(atoms::d_x_spin);
         cudaFree(atoms::d_spin);
         /*
         cu::atoms::x_spin_array.cu_real_array_t::~cu_real_array_t ();
         cu::atoms::y_spin_array.cu_real_array_t::~cu_real_array_t ();
         cu::atoms::z_spin_array.cu_real_array_t::~cu_real_array_t ();

         cu::atoms::x_coord_array.cu_real_array_t::~cu_real_array_t ();
         cu::atoms::y_coord_array.cu_real_array_t::~cu_real_array_t ();
         cu::atoms::z_coord_array.cu_real_array_t::~cu_real_array_t ();

         cu::atoms::type_array.cu_index_array_t::~cu_index_array_t ();
         cu::atoms::cell_array.cu_index_array_t::~cu_index_array_t ();
         cu::atoms::limits.cu_index_array_t::~cu_index_array_t ();
         cu::atoms::neighbours.cu_index_array_t::~cu_index_array_t ();

         cu::atoms::spin_norm_array.cu_real_array_t::~cu_real_array_t ();

         cu::cells::x_coord_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::y_coord_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::z_coord_array.cu_real_array_t::~cu_real_array_t ();

         cu::cells::x_mag_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::y_mag_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::z_mag_array.cu_real_array_t::~cu_real_array_t ();

         cu::cells::x_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::y_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::z_field_array.cu_real_array_t::~cu_real_array_t ();

         cu::cells::volume_array.cu_real_array_t::~cu_real_array_t ();
         cu::cells::num_atoms.cu_index_array_t::~cu_index_array_t ();

         cu::mp::materials.cu_material_array_t::~cu_material_array_t ();

         cu::x_total_spin_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::y_total_spin_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::z_total_spin_field_array.cu_real_array_t::~cu_real_array_t ();

         cu::x_total_external_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::y_total_external_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::z_total_external_field_array.cu_real_array_t::~cu_real_array_t ();

         cu::x_dipolar_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::y_dipolar_field_array.cu_real_array_t::~cu_real_array_t ();
         cu::z_dipolar_field_array.cu_real_array_t::~cu_real_array_t ();

         cu::llg::x_spin_buffer_array.cu_real_array_t::~cu_real_array_t ();
         cu::llg::y_spin_buffer_array.cu_real_array_t::~cu_real_array_t ();
         cu::llg::z_spin_buffer_array.cu_real_array_t::~cu_real_array_t ();

         cu::llg::dS_x_array.cu_real_array_t::~cu_real_array_t ();
         cu::llg::dS_y_array.cu_real_array_t::~cu_real_array_t ();
         cu::llg::dS_z_array.cu_real_array_t::~cu_real_array_t ();

         cu::llg::heun_parameters_device.thrust::device_vector<heun_parameters_t>::~device_vector ();

         cudaFree (d_rand_state);

         cu::stats::system_mask.cu_index_array_t::~cu_index_array_t ();
         cu::stats::system_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::system_mean_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::material_mask.cu_index_array_t::~cu_index_array_t ();
         cu::stats::material_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::material_mean_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::height_mask.cu_index_array_t::~cu_index_array_t ();
         cu::stats::height_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::height_mean_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::material_height_mask.cu_index_array_t::~cu_index_array_t ();
         cu::stats::material_height_magnetization.cu_real_array_t::~cu_real_array_t ();
         cu::stats::material_height_mean_magnetization.cu_real_array_t::~cu_real_array_t ();
         */

      }
   } /* internal */
#endif

} // end of namespace cuda
