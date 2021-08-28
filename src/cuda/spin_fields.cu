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
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "exchange_fields.hpp"
#include "data.hpp"
#include "internal.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

//------------------------------------------------------------------------------
// Host function to calculate spin fields using gpu kernel
//------------------------------------------------------------------------------
void update_spin_fields ()
{

   // Find the addresses in the device address space
   //int * d_materials = thrust::raw_pointer_cast(cu::atoms::type_array.data());
   //cu::material_parameters_t * d_material_params = thrust::raw_pointer_cast (cu::mp::materials.data());

   /*cu_real_t * d_x_spin = thrust::raw_pointer_cast(cu::atoms::x_spin_array.data());
   cu_real_t * d_y_spin = thrust::raw_pointer_cast(cu::atoms::y_spin_array.data());
   cu_real_t * d_z_spin = thrust::raw_pointer_cast(cu::atoms::z_spin_array.data());

   cu_real_t * d_x_spin_field = thrust::raw_pointer_cast(cu::x_total_spin_field_array.data());
   cu_real_t * d_y_spin_field = thrust::raw_pointer_cast(cu::y_total_spin_field_array.data());
   cu_real_t * d_z_spin_field = thrust::raw_pointer_cast(cu::z_total_spin_field_array.data());
   */
   // This exchange field zero out stuff, so they should come first
   cu::exchange::calculate_exchange_fields ();

   check_cuda_errors (__FILE__, __LINE__);

   // Call kernel to calculate non-exchange spin fields
   cu::update_non_exchange_spin_fields_kernel <<< cu::grid_size, cu::block_size >>> (
         cu::atoms::d_materials, cu::mp::d_material_params,
         cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
         cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
         ::atoms::num_atoms);

   check_cuda_errors (__FILE__, __LINE__);

}

//------------------------------------------------------------------------------
// Kernel function to calculate external fields
//------------------------------------------------------------------------------
__global__ void update_non_exchange_spin_fields_kernel (
      int * material,
      cu::material_parameters_t * material_params,
      cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
      cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
      int n_atoms
      )
{
   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x){

      int mid = material[i];
      // Load parameters from memory
      cu::material_parameters_t material = material_params[mid];

      // Initialise register to hold total spin field
      cu_real_t field_x = 0.0;
      cu_real_t field_y = 0.0;
      cu_real_t field_z = 0.0;

      cu_real_t sx = x_spin[i];
      cu_real_t sy = y_spin[i];
      cu_real_t sz = z_spin[i];

      cu_real_t ex = material.anisotropy_unit_x;
      cu_real_t ey = material.anisotropy_unit_y;
      cu_real_t ez = material.anisotropy_unit_z;

      cu_real_t sdote = sx * ex + sy * ey + sz * ez;
      cu_real_t sdote3 = sdote * sdote * sdote;
      cu_real_t sdote5 = sdote3 * sdote * sdote;

      /*
       * Spherical harmonics
       */

      cu_real_t scale = 0.6666666666666667;

      cu_real_t k2 = material.sh2;
      cu_real_t k4 = material.sh4;
      cu_real_t k6 = material.sh6;

      cu_real_t ek2 = k2 * 3.0 * sdote;
      cu_real_t ek4 = k4 * 0.125 * (140.0 * sdote3 - 60.0 *sdote);
      cu_real_t ek6 = k6 * 0.0625 * (1386.0 * sdote5 - 1260.0 * sdote3 + 210.0 * sdote);

      field_x += scale * ex * (ek2 + ek4 + ek6);
      field_y += scale * ey * (ek2 + ek4 + ek6);
      field_z += scale * ez * (ek2 + ek4 + ek6);

      /*
       * Lattice anisotropy
       */

      /*
       * TODO: add the temperature dependence
       */

      /*
       * TODO: communicate every timestep
       */

      //cu_real_t k_latt = 2.0 * material.k_latt;
      //field_x -= k_latt * ex * sdote;
      //field_y -= k_latt * ey * sdote;
      //field_z -= k_latt * ez * sdote;

      /*
       * TODO: Surface anisotropy?
       */

      /*
       * TODO: Lagrange multipliers?
       */

      /*
       * Write back to main memory
       */

      x_sp_field[i] += field_x;
      y_sp_field[i] += field_y;
      z_sp_field[i] += field_z;

   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
