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
#include "sim.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "exchange_fields.hpp"
#include "data.hpp"
#include "dipole.hpp"
#include "internal.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

void update_dipolar_fields ()
{

   // check if dipole calculation is enabled
   if(sim::hamiltonian_simulation_flags[4]!=1) return;

   // check for previous demag update at same time (avoids recalculation in Heun scheme)
   if (::sim::time == ::dipole::update_time) return;

   // if remainder of time/rate != 0 return
   if (::sim::time % ::dipole::update_rate != 0) return;

   // save last time of demag update
   ::dipole::update_time = ::sim::time;

   update_cell_magnetizations ();

   check_cuda_errors (__FILE__, __LINE__);

   /*
    * Figure out addresses in device memory space
    */

   cu_real_t * d_x_mag = thrust::raw_pointer_cast(
         cu::cells::x_mag_array.data());
   cu_real_t * d_y_mag = thrust::raw_pointer_cast(
         cu::cells::y_mag_array.data());
   cu_real_t * d_z_mag = thrust::raw_pointer_cast(
         cu::cells::z_mag_array.data());

   cu_real_t * d_x_coord = thrust::raw_pointer_cast(
         cu::cells::x_coord_array.data());
   cu_real_t * d_y_coord = thrust::raw_pointer_cast(
         cu::cells::y_coord_array.data());
   cu_real_t * d_z_coord = thrust::raw_pointer_cast(
         cu::cells::z_coord_array.data());

   cu_real_t * d_volume = thrust::raw_pointer_cast(
         cu::cells::volume_array.data());

   cu_real_t * d_x_cell_field = thrust::raw_pointer_cast(
         cu::cells::x_field_array.data());
   cu_real_t * d_y_cell_field = thrust::raw_pointer_cast(
         cu::cells::y_field_array.data());
   cu_real_t * d_z_cell_field = thrust::raw_pointer_cast(
         cu::cells::z_field_array.data());

   /*
    * Update cell dipolar fields
    */

   update_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
         d_x_mag, d_y_mag, d_z_mag,
         d_x_coord, d_y_coord, d_z_coord,
         d_volume,
         d_x_cell_field, d_y_cell_field, d_z_cell_field,
         ::cells::num_cells
         );

   check_cuda_errors (__FILE__, __LINE__);

   /*
    * Update atomistic dipolar fields
    */

   int * d_cells =
      thrust::raw_pointer_cast(cu::atoms::cell_array.data());

   cu_real_t * d_x_atom_field = thrust::raw_pointer_cast(
         cu::x_dipolar_field_array.data());
   cu_real_t * d_y_atom_field = thrust::raw_pointer_cast(
         cu::y_dipolar_field_array.data());
   cu_real_t * d_z_atom_field = thrust::raw_pointer_cast(
         cu::z_dipolar_field_array.data());

   update_atomistic_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
         d_x_cell_field, d_y_cell_field, d_z_cell_field,
         d_x_atom_field, d_y_atom_field, d_z_atom_field,
         d_cells,
         ::atoms::num_atoms
         );

   check_cuda_errors (__FILE__, __LINE__);
}

void update_cell_magnetizations ()
{
   cu_real_t * d_x_spin = thrust::raw_pointer_cast(
         cu::atoms::x_spin_array.data());
   cu_real_t * d_y_spin = thrust::raw_pointer_cast(
         cu::atoms::y_spin_array.data());
   cu_real_t * d_z_spin = thrust::raw_pointer_cast(
         cu::atoms::z_spin_array.data());

   int * d_materials =
      thrust::raw_pointer_cast(cu::atoms::type_array.data());

   int * d_cells =
      thrust::raw_pointer_cast(cu::atoms::cell_array.data());

   cu::material_parameters_t * d_material_params =
      thrust::raw_pointer_cast (cu::mp::materials.data());

   cu_real_t * d_x_mag = thrust::raw_pointer_cast(
         cu::cells::x_mag_array.data());
   cu_real_t * d_y_mag = thrust::raw_pointer_cast(
         cu::cells::y_mag_array.data());
   cu_real_t * d_z_mag = thrust::raw_pointer_cast(
         cu::cells::z_mag_array.data());

   /*
    * Update cell magnetizations
    */

   thrust::fill(
         cu::cells::x_mag_array.begin(),
         cu::cells::x_mag_array.end(),
         0.0);
   thrust::fill(
         cu::cells::y_mag_array.begin(),
         cu::cells::y_mag_array.end(),
         0.0);
   thrust::fill(
         cu::cells::z_mag_array.begin(),
         cu::cells::z_mag_array.end(),
         0.0);

   update_cell_magnetization <<< cu::grid_size, cu::block_size >>> (
         d_x_spin, d_y_spin, d_z_spin,
         d_materials, d_cells,
         d_material_params,
         d_x_mag, d_y_mag, d_z_mag,
         ::atoms::num_atoms
         );

   check_cuda_errors (__FILE__, __LINE__);
}

__global__ void update_cell_magnetization (
      cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
      int * material, int * cell,
      cu::material_parameters_t * material_params,
      cu_real_t * x_mag, cu_real_t * y_mag, cu_real_t * z_mag,
      int n_atoms
      )
{
   /*
    * TODO: This is an supremely naïve implementation
    *       the number of cells can be as big as the number of atoms
    *       so might as well leave it like this
    */

   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x)
   {
      int mid = material[i];
      int cid = cell[i];
      cu_real_t mu_s = material_params[mid].mu_s_si;
      cu::atomicAdd(&x_mag[cid], x_spin[i] * mu_s);
      cu::atomicAdd(&y_mag[cid], y_spin[i] * mu_s);
      cu::atomicAdd(&z_mag[cid], z_spin[i] * mu_s);
   }
}

__global__ void update_dipolar_fields (
      cu_real_t * x_mag, cu_real_t * y_mag, cu_real_t * z_mag,
      cu_real_t * x_coord, cu_real_t * y_coord, cu_real_t * z_coord,
      cu_real_t * volume,
      cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
      int n_cells
      )
{
   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_cells;
         i += blockDim.x * gridDim.x)
   {
      cu_real_t mx = x_mag[i];
      cu_real_t my = y_mag[i];
      cu_real_t mz = z_mag[i];
      cu_real_t cx = x_coord[i];
      cu_real_t cy = y_coord[i];
      cu_real_t cz = z_coord[i];
      /*
      * Inverse volume from the number of atoms in macro-cell
      */
      cu_real_t vol_prefac = - 4.0 * M_PI / (3.0 * volume[i]);
      cu_real_t prefactor = 1.0e+23; // 1e-7/1e30

      cu_real_t field_x = vol_prefac * mx;
      cu_real_t field_y = vol_prefac * my;
      cu_real_t field_z = vol_prefac * mz;

      for (int j = 0; j < n_cells; j++)
      {
         if (i == j) continue;
         cu_real_t omx = x_mag[i];
         cu_real_t omy = y_mag[i];
         cu_real_t omz = z_mag[i];

         cu_real_t dx = x_coord[j] - cx;
         cu_real_t dy = y_coord[j] - cy;
         cu_real_t dz = z_coord[j] - cz;

         cu_real_t drij = rsqrt(dx * dx + dy * dy + dz * dz);
         cu_real_t drij3 = drij * drij * drij;

         cu_real_t sdote = (
               omx * dx * drij +
               omy * dy * drij +
               omz * dz * drij);

         field_x += (3.0 * sdote * dx * drij - omx) * drij3;
         field_y += (3.0 * sdote * dy * drij - omy) * drij3;
         field_z += (3.0 * sdote * dz * drij - omz) * drij3;
      }

      x_dip_field[i] = prefactor * field_x;
      y_dip_field[i] = prefactor * field_y;
      z_dip_field[i] = prefactor * field_z;
   }
}

__global__ void update_atomistic_dipolar_fields (
      cu_real_t * x_cell_field, cu_real_t * y_cell_field, cu_real_t * z_cell_field,
      cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
      int * cell, int n_atoms
      )
{
   /*
    * TODO: Warning extremely naïve data access pattern
    */
   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x)
   {
      int cid = cell[i];
      x_dip_field[i] = x_cell_field[cid];
      y_dip_field[i] = y_cell_field[cid];
      z_dip_field[i] = z_cell_field[cid];
   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
