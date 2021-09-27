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
