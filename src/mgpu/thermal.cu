//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) S Jenkins 2026. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cstdint>

// CUDA headers
#include <curand_kernel.h>

// Vampire headers
#include "hamr.hpp"
#include "mgpu.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// CUDA module headers
#include "../cuda/cuda_utils.hpp"
#include "../cuda/data.hpp"
#include "../cuda/internal.hpp"
#include "../cuda/typedefs.hpp"

// mgpu module headers
#include "internal.hpp"

namespace cu = ::vcuda::internal;

//------------------------------------------------------------------------------
// Thermal and HAMR fields as in src/cuda/thermal_fields.cu and hamr_fields.cu,
// except that the random numbers come from the Philox generator seeded by the
// integration seed, atom key and step number, so each atom has the same noise
// on any number of GPUs.
//------------------------------------------------------------------------------
namespace mgpu{

namespace internal{

   uint64_t* d_atom_key = nullptr;
   uint64_t noise_step = 0; // number of the next thermal field update

   //---------------------------------------------------------------------------
   // Three independent standard normal random numbers for one atom and step,
   // from one Philox4x32-10 draw and cuRAND's Box-Muller transform (single
   // precision, as the CUDA module)
   //---------------------------------------------------------------------------
   __device__ float3 thermal_noise(const uint64_t key, const uint64_t step, const uint2 seed){
      const uint4 counter = make_uint4(uint32_t(key), uint32_t(key >> 32), uint32_t(step), uint32_t(step >> 32));
      const uint4 r = curand_Philox4x32_10(counter, seed);
      const float2 a = _curand_box_muller(r.x, r.y);
      const float2 b = _curand_box_muller(r.z, r.w);
      return make_float3(a.x, a.y, b.x);
   }

   //---------------------------------------------------------------------------
   // Thermal field at a uniform temperature (cf. apply_global_temperature_kernel)
   //---------------------------------------------------------------------------
   __global__ void global_temperature_kernel(
         cu::cu_real_t * x_field_array, cu::cu_real_t * y_field_array, cu::cu_real_t * z_field_array,
         const cu::cu_real_t temperature,
         const uint64_t * atom_key, const uint64_t step, const uint2 seed,
         const cu::material_parameters_t * material_params,
         const int * material,
         const int n_atoms,
         const int * active){

      for(int k = blockIdx.x * blockDim.x + threadIdx.x; k < n_atoms; k += blockDim.x * gridDim.x){

         const int i = active ? active[k] : k; // all atoms, or the listed ones (HAMR active window)

         const cu::material_parameters_t mat = material_params[material[i]];
         if(mat.skip) continue; // atoms that are not integrated need no field (their noise depends on no other atom)
         const cu::cu_real_t rsigma = cu::thermal_field_sigma(mat, temperature);

         const float3 n = thermal_noise(atom_key[i], step, seed);
         x_field_array[i] = rsigma * n.x;
         y_field_array[i] = rsigma * n.y;
         z_field_array[i] = rsigma * n.z;

      }

   }

   //---------------------------------------------------------------------------
   // Thermal field of the HAMR laser profile and the head field inside the
   // field box (cf. apply_local_temperature_and_field_kernel)
   //---------------------------------------------------------------------------
   __global__ void local_temperature_and_field_kernel(
         const cu::cu_real_t * atoms_coord_x, const cu::cu_real_t * atoms_coord_y,
         cu::cu_real_t * x_field_array, cu::cu_real_t * y_field_array, cu::cu_real_t * z_field_array,
         const cu::cu_real_t Tmin, const cu::cu_real_t DeltaT,
         const cu::cu_real_t one_over_denx, const cu::cu_real_t one_over_deny,
         const cu::cu_real_t px, const cu::cu_real_t py,
         const cu::cu_real_t Hx_app, const cu::cu_real_t Hy_app, const cu::cu_real_t Hz_app,
         const cu::cu_real_t Hloc_min_x, const cu::cu_real_t Hloc_max_x,
         const cu::cu_real_t Hloc_min_y, const cu::cu_real_t Hloc_max_y,
         const uint64_t * atom_key, const uint64_t step, const uint2 seed,
         const cu::material_parameters_t * material_params,
         const int * material,
         const int n_atoms,
         const int * active){

      for(int k = blockIdx.x * blockDim.x + threadIdx.x; k < n_atoms; k += blockDim.x * gridDim.x){

         const int i = active ? active[k] : k; // all atoms, or the listed ones (HAMR active window)

         const cu::material_parameters_t mat = material_params[material[i]];
         if(mat.skip) continue; // atoms that are not integrated need no field (their noise depends on no other atom)

         // temperature of the Gaussian laser profile and thermal field with the material's temperature rescaling
         const cu::cu_real_t cx = atoms_coord_x[i];
         const cu::cu_real_t cy = atoms_coord_y[i];
         const cu::cu_real_t temp = cu::laser_temperature(cx, cy, px, py, Tmin, DeltaT, one_over_denx, one_over_deny);
         const cu::cu_real_t rsigma = cu::thermal_field_sigma(mat, temp);

         const float3 n = thermal_noise(atom_key[i], step, seed);
         cu::cu_real_t field_x = rsigma * n.x;
         cu::cu_real_t field_y = rsigma * n.y;
         cu::cu_real_t field_z = rsigma * n.z;

         // add the head field if the atom is within the field box
         if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
            field_x += Hx_app;
            field_y += Hy_app;
            field_z += Hz_app;
         }

         x_field_array[i] = field_x;
         y_field_array[i] = field_y;
         z_field_array[i] = field_z;

      }

   }

   //---------------------------------------------------------------------------
   // Function to copy the atom keys to the device and to start counting thermal
   // field updates at the current time step (so that a simulation continued
   // from a checkpoint does not repeat the noise of its first steps)
   //---------------------------------------------------------------------------
   void initialize_noise(){
      cudaMalloc((void**)&d_atom_key, atom_key.size() * sizeof(uint64_t));
      cudaMemcpy(d_atom_key, atom_key.data(), atom_key.size() * sizeof(uint64_t), cudaMemcpyHostToDevice);
      noise_step = sim::time;
   }

   void free_noise(){
      cudaFree(d_atom_key);
      d_atom_key = nullptr;
   }

   //---------------------------------------------------------------------------
   // Function to calculate the external fields of one time step: as the CUDA
   // module's update_external_fields, with the thermal and HAMR fields above
   //---------------------------------------------------------------------------
   void update_external_fields(){

      // halo atoms are not integrated, so they need no thermal and HAMR fields
      const int num_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      const uint2 seed = make_uint2(uint32_t(::mtrandom::integration_seed), 0u);
      const uint64_t step = noise_step++;

      // HAMR simulation: laser profile and head field, or the global temperature with the laser off
      if(::program::program == ::program::program_t::hamr_simulation){
         if(::hamr::head_laser_on){
            const cu::cu_real_t Tmin = sim::Tmin;
            const cu::cu_real_t DeltaT = sim::Tmax - sim::Tmin;
            const cu::cu_real_t laser_sigma_x2 = cu::hamr::d_laser_sigma_x * cu::hamr::d_laser_sigma_x;
            const cu::cu_real_t laser_sigma_y2 = cu::hamr::d_laser_sigma_y * cu::hamr::d_laser_sigma_y;
            const cu::cu_real_t one_over_denx = 1.0/(2.0 * laser_sigma_x2);
            const cu::cu_real_t one_over_deny = 1.0/(2.0 * laser_sigma_y2);
            const cu::cu_real_t Hx_app = sim::H_vec[0]*sim::H_applied;
            const cu::cu_real_t Hy_app = sim::H_vec[1]*sim::H_applied;
            const cu::cu_real_t Hz_app = sim::H_vec[2]*sim::H_applied;

            // laser and field box centred on the current head position
            cu::hamr::d_head_position_x = ::hamr::get_head_position_x();
            cu::hamr::d_head_position_y = ::hamr::get_head_position_y();
            const cu::cu_real_t px = cu::hamr::d_head_position_x;
            const cu::cu_real_t py = cu::hamr::d_head_position_y;
            const cu::cu_real_t Hloc_min_x = px - 0.5*cu::hamr::d_H_bounds_x - cu::hamr::d_NPS;
            const cu::cu_real_t Hloc_max_x = px + 0.5*cu::hamr::d_H_bounds_x - cu::hamr::d_NPS;
            const cu::cu_real_t Hloc_min_y = py - 0.5*cu::hamr::d_H_bounds_y;
            const cu::cu_real_t Hloc_max_y = py + 0.5*cu::hamr::d_H_bounds_y;

            local_temperature_and_field_kernel <<< cu::grid_size, cu::block_size >>> (
               cu::atoms::d_x_coord, cu::atoms::d_y_coord,
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
               Tmin, DeltaT, one_over_denx, one_over_deny, px, py,
               Hx_app, Hy_app, Hz_app, Hloc_min_x, Hloc_max_x, Hloc_min_y, Hloc_max_y,
               d_atom_key, step, seed, cu::mp::d_material_params, cu::atoms::d_materials,
               cu::num_integrated_atoms(num_atoms), cu::d_active_atoms);
         }
         else{
            global_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
               sim::temperature, d_atom_key, step, seed, cu::mp::d_material_params, cu::atoms::d_materials,
               cu::num_integrated_atoms(num_atoms), cu::d_active_atoms);
         }
      }
      // otherwise the global thermal field and the applied fields
      else{
         global_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
            cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
            sim::temperature, d_atom_key, step, seed, cu::mp::d_material_params, cu::atoms::d_materials, num_atoms, nullptr);
         cu::update_applied_fields();
      }
      check_cuda_errors(__FILE__, __LINE__);

      // oscillating field once enabled by the FMR program
      if(sim::enable_fmr) cu::update_fmr_fields();

   }

} // end of internal namespace

} // end of mgpu namespace
