/**
 * @brief this file provides definitions for the off-topic internal function
 *        definitions.
 */

#include "data.hpp"
#include "internal.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda
{
#ifdef CUDA
   namespace internal
   {
      bool __finalize ()
      {
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

         cu::cells::x_coord_array.RealArray::~device_vector ();
         cu::cells::y_coord_array.RealArray::~device_vector ();
         cu::cells::z_coord_array.RealArray::~device_vector ();
         cu::cells::x_mag_array.RealArray::~device_vector ();
         cu::cells::y_mag_array.RealArray::~device_vector ();
         cu::cells::z_mag_array.RealArray::~device_vector ();
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

         return true;
      }
   } /* internal */
#endif
} /* vcuda */
