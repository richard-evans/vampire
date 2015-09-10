//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "../../hdr/cuda.hpp"

// Local cuda headers
#include "data.hpp"
#include "llg_heun.hpp"
#include "internal.hpp"

#include "exchange_fields.hpp"
#include "statistics.hpp"

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

#endif

      return;
   }

#ifdef CUDA
   namespace internal
   {
      void __finalize ()
      {

         // De-allocate the exchange fields
         cu::exchange::finalise_exchange();

         cu::atoms::x_spin_array.RealArray::~device_vector ();
         cu::atoms::y_spin_array.RealArray::~device_vector ();
         cu::atoms::z_spin_array.RealArray::~device_vector ();
         cu::atoms::x_coord_array.RealArray::~device_vector ();
         cu::atoms::y_coord_array.RealArray::~device_vector ();
         cu::atoms::z_coord_array.RealArray::~device_vector ();
         cu::atoms::type_array.IndexArray::~device_vector ();
         cu::atoms::cell_array.IndexArray::~device_vector ();
         cu::atoms::limits.IndexArray::~device_vector ();
         cu::atoms::neighbours.IndexArray::~device_vector ();
         cu::atoms::spin_norm_array.RealArray::~device_vector ();

         exchange::Jxx_vals_d.RealArray::~device_vector ();
         exchange::Jyy_vals_d.RealArray::~device_vector ();
         exchange::Jzz_vals_d.RealArray::~device_vector ();

         cu::cells::x_coord_array.RealArray::~device_vector ();
         cu::cells::y_coord_array.RealArray::~device_vector ();
         cu::cells::z_coord_array.RealArray::~device_vector ();
         cu::cells::x_mag_array.RealArray::~device_vector ();
         cu::cells::y_mag_array.RealArray::~device_vector ();
         cu::cells::z_mag_array.RealArray::~device_vector ();
         cu::cells::x_field_array.RealArray::~device_vector ();
         cu::cells::y_field_array.RealArray::~device_vector ();
         cu::cells::z_field_array.RealArray::~device_vector ();
         cu::cells::volume_array.RealArray::~device_vector ();
         cu::cells::num_atoms.IndexArray::~device_vector ();

         cu::mp::materials.MaterialParametersArray::~device_vector ();

         cu::x_total_spin_field_array.RealArray::~device_vector ();
         cu::y_total_spin_field_array.RealArray::~device_vector ();
         cu::z_total_spin_field_array.RealArray::~device_vector ();
         cu::x_total_external_field_array.RealArray::~device_vector ();
         cu::y_total_external_field_array.RealArray::~device_vector ();
         cu::z_total_external_field_array.RealArray::~device_vector ();
         cu::x_dipolar_field_array.RealArray::~device_vector ();
         cu::y_dipolar_field_array.RealArray::~device_vector ();
         cu::z_dipolar_field_array.RealArray::~device_vector ();

         cu::llg::x_spin_buffer_array.RealArray::~device_vector ();
         cu::llg::y_spin_buffer_array.RealArray::~device_vector ();
         cu::llg::z_spin_buffer_array.RealArray::~device_vector ();

         cu::llg::heun_parameters.HeunParametersArray::~device_vector ();

         cudaFree (d_rand_state);

         cu::stats::system_mask.IndexArray::~device_vector ();
         cu::stats::system_magnetization.RealArray::~device_vector ();
         cu::stats::system_mean_magnetization.RealArray::~device_vector ();
         cu::stats::material_mask.IndexArray::~device_vector ();
         cu::stats::material_magnetization.RealArray::~device_vector ();
         cu::stats::material_mean_magnetization.RealArray::~device_vector ();
         cu::stats::height_mask.IndexArray::~device_vector ();
         cu::stats::height_magnetization.RealArray::~device_vector ();
         cu::stats::height_mean_magnetization.RealArray::~device_vector ();
         cu::stats::material_height_mask.IndexArray::~device_vector ();
         cu::stats::material_height_magnetization.RealArray::~device_vector ();
         cu::stats::material_height_mean_magnetization.RealArray::~device_vector ();

      }
   } /* internal */
#endif

} // end of namespace cuda