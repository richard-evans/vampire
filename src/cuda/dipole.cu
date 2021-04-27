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
    std::cout << "update_dipolar_field() called" << std::endl;

   // check if dipole calculation is enabled
   if(!::dipole::activated) return;

   // check for previous demag update at same time (avoids recalculation in Heun scheme)
   if (::sim::time == ::dipole::update_time) return;

   // if remainder of time/rate != 0 return
   if (::sim::time % ::dipole::update_rate != 0) return;

   // save last time of demag update
   ::dipole::update_time = ::sim::time;

   update_cell_magnetizations ();

   check_cuda_errors (__FILE__, __LINE__);

   int num_local_cells = ::dipole::get_tot_num_local_cells();
   int num_cells = ::dipole::get_tot_num_cells();

   int *d_cell_id_array;
   cudaMalloc((void**)&d_cell_id_array, ::cells::cell_id_array.size() * sizeof(int));
   cudaMemcpy(d_cell_id_array, ::cells::cell_id_array.data(), ::cells::cell_id_array.size() * sizeof(int), cudaMemcpyHostToDevice);

   std::vector<int> num_atoms_in_cell = ::dipole::get_num_atoms_in_cell_array();
   int *d_num_atoms_in_cell;
   cudaMalloc((void**)&d_num_atoms_in_cell, num_atoms_in_cell.size() * sizeof(int));
   cudaMemcpy(d_num_atoms_in_cell, num_atoms_in_cell.data(), num_atoms_in_cell.size() * sizeof(int), cudaMemcpyHostToDevice);


   /*
    * Figure out addresses in device memory space
    */

   /*cu_real_t * d_x_mag = thrust::raw_pointer_cast(
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
*/
   /*
    * Update cell dipolar fields
    */

   update_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
         cu::cells::d_x_mag, cu::cells::d_y_mag, cu::cells::d_z_mag,
         cu::cells::d_x_coord, cu::cells::d_y_coord, cu::cells::d_z_coord,
         cu::cells::d_volume,
         cu::cells::d_x_cell_field, cu::cells::d_y_cell_field, cu::cells::d_z_cell_field,
         cu::cells::d_x_cell_mu0H_field, cu::cells::d_y_cell_mu0H_field, cu::cells::d_z_cell_mu0H_field,
         cu::cells::d_tensor_xx, cu::cells::d_tensor_xy, cu::cells::d_tensor_xz,
         cu::cells::d_tensor_yy, cu::cells::d_tensor_yz, cu::cells::d_tensor_zz,
         d_cell_id_array,
         d_num_atoms_in_cell,
         num_local_cells,
         num_cells
//         ::cells::num_cells
         );

   check_cuda_errors (__FILE__, __LINE__);

   /*
    * Update atomistic dipolar fields
    */

    /*
   int * d_cells =
      thrust::raw_pointer_cast(cu::atoms::cell_array.data());

   cu_real_t * d_x_atom_field = thrust::raw_pointer_cast(
         cu::x_dipolar_field_array.data());
   cu_real_t * d_y_atom_field = thrust::raw_pointer_cast(
         cu::y_dipolar_field_array.data());
   cu_real_t * d_z_atom_field = thrust::raw_pointer_cast(
         cu::z_dipolar_field_array.data());
      */
   update_atomistic_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
         cu::cells::d_x_cell_field, cu::cells::d_y_cell_field, cu::cells::d_z_cell_field,
         cu::d_x_dip_field, cu::d_y_dip_field, cu::d_z_dip_field,
         cu::atoms::d_cells,
         ::atoms::num_atoms
         );

   check_cuda_errors (__FILE__, __LINE__);
   
   // Transfer cells dipolar fields from gpu to cpu
   vcuda::transfer_dipole_cells_fields_from_gpu_to_cpu();
}

void update_cell_magnetizations ()
{
   /*
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
   */
   /*
    * Update cell magnetizations
    */

    cudaMemset(cu::cells::d_x_mag, 0, ::cells::num_cells * sizeof(cu_real_t));
    cudaMemset(cu::cells::d_y_mag, 0, ::cells::num_cells * sizeof(cu_real_t));
    cudaMemset(cu::cells::d_z_mag, 0, ::cells::num_cells * sizeof(cu_real_t));

/*
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
*/
   update_cell_magnetization <<< cu::grid_size, cu::block_size >>> (
         cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
         cu::atoms::d_materials, cu::atoms::d_cells,
         cu::mp::d_material_params,
         cu::cells::d_x_mag, cu::cells::d_y_mag, cu::cells::d_z_mag,
         ::atoms::num_atoms
         );

   check_cuda_errors (__FILE__, __LINE__);
}

// Some data is read only (material, cell, material_params)
// and access in non-coalesced manner - should benefit from const __restrict__
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
      // What are the ranges and general layout of mid and cid
      // Could stuff be moved to shared memory and
      // warp reduction used to avoid atomicAdd?
      int mid = material[i];
      int cid = cell[i];
      cu_real_t mu_s = material_params[mid].mu_s_si;
      cu::atomicAdd(&x_mag[cid], x_spin[i] * mu_s);
      cu::atomicAdd(&y_mag[cid], y_spin[i] * mu_s);
      cu::atomicAdd(&z_mag[cid], z_spin[i] * mu_s);
   }
}

// Same as above - some data is read only, although looks like it's
// at least accessed in contoguous manner
__global__ void update_dipolar_fields (
      cu_real_t * x_mag, cu_real_t * y_mag, cu_real_t * z_mag,
      cu_real_t * x_coord, cu_real_t * y_coord, cu_real_t * z_coord,
      cu_real_t * volume,
      cu_real_t * x_cell_field, cu_real_t * y_cell_field, cu_real_t * z_cell_field,
      cu_real_t * x_cell_mu0H_field, cu_real_t * y_cell_mu0H_field, cu_real_t * z_cell_mu0H_field,
      cu_real_t * d_tensor_xx, cu_real_t * d_tensor_xy, cu_real_t * d_tensor_xz,
      cu_real_t * d_tensor_yy, cu_real_t * d_tensor_yz, cu_real_t * d_tensor_zz,
      int * d_cell_id_array,
      int * d_num_atoms_in_cell,
      int n_local_cells, int n_cells
      )
{

   cu_real_t imuB = 1.0/9.27400915e-24;
   cu_real_t prefactor = 9.27400915e-01;     // prefactor = mu_B * (mu_0/(4*pi) /1e-30)

   // Define counter for 1D dipole tensor
   int i_1Dindex = 0;

   for ( int lc = blockIdx.x * blockDim.x + threadIdx.x;
         lc < n_local_cells;
         lc += blockDim.x * gridDim.x)
   {

      int i = d_cell_id_array[lc]; //::cells::cell_id_array[lc];

      if(d_num_atoms_in_cell[i]>0){

         /*
         cu_real_t cx = x_coord[i];
         cu_real_t cy = y_coord[i];
         cu_real_t cz = z_coord[i];
         */

         /*
         * Inverse volume from the number of atoms in macro-cell
         */
         // DEFINES ? 4.0 * PI / 3 is constant
         // Should be optimised out anyway?
         /*
         cu_real_t vol_prefac = -4.0 * M_PI / (3.0 * volume[i]);
         cu_real_t prefactor = 1.0e+23; // 1e-7/1e30
         */

         cu_real_t self_demag = 8.0 * M_PI / (3.0 * volume[i]);

         // Normalise cells magnetisation
         cu_real_t mx_i = x_mag[i] * imuB;
         cu_real_t my_i = y_mag[i] * imuB;
         cu_real_t mz_i = z_mag[i] * imuB;

         // Add self-demagnetisation term
         cu_real_t field_x = self_demag * mx_i * 0.0;
         cu_real_t field_y = self_demag * my_i * 0.0;
         cu_real_t field_z = self_demag * mz_i * 0.0;

         // Add self-demagnetisation term
         cu_real_t mu0Hd_field_x = -0.5 * self_demag * mx_i;
         cu_real_t mu0Hd_field_y = -0.5 * self_demag * my_i;
         cu_real_t mu0Hd_field_z = -0.5 * self_demag * mz_i;


         for ( int j = blockIdx.x * blockDim.x + threadIdx.x;
              j < n_cells;
              j += blockDim.x * gridDim.x)
         {
            if(d_num_atoms_in_cell[j]>0){

               cu_real_t mx_j = x_mag[j] * imuB;
               cu_real_t my_j = y_mag[j] * imuB;
               cu_real_t mz_j = z_mag[j] * imuB;

               field_x += (mx_j * d_tensor_xx[i_1Dindex] + my_j * d_tensor_xy[i_1Dindex] + mz_j * d_tensor_xz[i_1Dindex]);
               field_y += (mx_j * d_tensor_xy[i_1Dindex] + my_j * d_tensor_yy[i_1Dindex] + mz_j * d_tensor_yz[i_1Dindex]);
               field_z += (mx_j * d_tensor_xz[i_1Dindex] + my_j * d_tensor_yz[i_1Dindex] + mz_j * d_tensor_zz[i_1Dindex]);

               mu0Hd_field_x += (mx_j * d_tensor_xx[i_1Dindex] + my_j * d_tensor_xy[i_1Dindex] + mz_j * d_tensor_xz[i_1Dindex]);
               mu0Hd_field_y += (mx_j * d_tensor_xy[i_1Dindex] + my_j * d_tensor_yy[i_1Dindex] + mz_j * d_tensor_yz[i_1Dindex]);
               mu0Hd_field_z += (mx_j * d_tensor_xz[i_1Dindex] + my_j * d_tensor_yz[i_1Dindex] + mz_j * d_tensor_zz[i_1Dindex]);

               /*
               // Make use of float3?
               // Store in AoS instead of separate arrays?
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
               */
            }
            else{ // Increase counter if cell j is empty
               i_1Dindex++;
            } // end if cell i is not empty
         } // end for loop over n_cells

         // Same AoS argument as above?
         x_cell_field[i] = prefactor * field_x;
         y_cell_field[i] = prefactor * field_y;
         z_cell_field[i] = prefactor * field_z;

         x_cell_mu0H_field[i] = prefactor * mu0Hd_field_x;
         y_cell_mu0H_field[i] = prefactor * mu0Hd_field_y;
         z_cell_mu0H_field[i] = prefactor * mu0Hd_field_z;
      }
      else{ // Increase counter of n_cells if cell i is empty to go to next i
         i_1Dindex += n_cells;
      } // end if cell j is not empty
   } // end for loop over n_local_cells
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
      // I bet this is non-coalescing?
      // Can it somehow be made coalescing with shared memory?
      int cid = cell[i];
      x_dip_field[i] = x_cell_field[cid];
      y_dip_field[i] = y_cell_field[cid];
      z_dip_field[i] = z_cell_field[cid];
   }
}

} // end of internal namespace

void update_dipolar_fields (){
//    std::cout << "Updating dipolar filed with GPU" << std::endl;
    cu::update_dipolar_fields();
}

} // end of vcuda namespace

#endif
