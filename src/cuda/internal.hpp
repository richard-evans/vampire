#ifndef CUDA_INTERNAL_H_
#define CUDA_INTERNAL_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// cuda implementation. These functions should not be accessed outside of the
// cuda code.
//---------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <vector>

// Include type definitions for cuda code
#include "typedefs.hpp"

/*
 * requesting data strcutures from the main program
 */

#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"

#include "cuda_timer.h"

namespace vcuda {

#ifdef CUDA

namespace internal {

extern timer cuda_timer;

/*
 * Thread launch parameters
 */

extern int block_size;
extern int grid_size;

/*
 * Internal data structures
 */

/*
 * Initlialization functions
 */
bool __initialize_atoms();
bool __initialize_fields();
bool __initialize_neel();
bool __initialize_cells();
bool __initialize_dipole();
bool __initialize_hamr();
bool __initialize_materials();
bool __initialize_materials_rotational();
bool __initialize_topology();
bool __initialize_curand();
bool __initialize_stats();

/*
 * Clean up function
 */
void __finalize();

/*
 * Field updates
 */

void update_spin_fields(const int num_atoms, const int first_atom = 0); // atoms first_atom to num_atoms-1
void update_external_fields();
void update_dipolar_fields();
void update_cell_magnetizations();
void update_hamr_field();
// atoms integrated with the HAMR active window (hamr:active-window), in increasing order; nullptr: all atoms
extern int *d_active_atoms;
extern std::vector<int> h_active_atoms; // host copy
void update_active_atoms(); // copies the window's atoms to d_active_atoms when it changes
extern std::vector<int> active_ranges; // the window's atoms as [first, last+1) pairs
extern int active_window_version;      // incremented when the window changes
// number of the atoms 0 to num_atoms-1 the kernels loop over: all, or those in the active window; kernels
// for the atoms first to last-1 loop over positions num_integrated_atoms(first) to num_integrated_atoms(last)-1
inline int num_integrated_atoms(const int num_atoms){
   return d_active_atoms ? std::lower_bound(h_active_atoms.begin(), h_active_atoms.end(), num_atoms) - h_active_atoms.begin() : num_atoms;
}
void update_global_thermal_field();
void update_applied_fields();
void update_fmr_fields();
void update_neel_fields(const int num_atoms, const int first_atom = 0);
void update_rotational_fields(const int num_atoms, const int first_atom = 0);

/*
 * Shared device functions
 */

inline __device__ double atomicAdd(double *address, double value) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull;
  unsigned long long int assumed;
  do {
    assumed = old;
    old =
        atomicCAS(address_as_ull, assumed,
                  __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

inline __device__ float atomicAdd(float *address, float value) {
  return ::atomicAdd(address, value);
}

// width of the thermal field of a material at a temperature (with the material's temperature
// rescaling); shared by the thermal field kernels of the CUDA and multi-GPU versions
inline __device__ cu_real_t thermal_field_sigma(const material_parameters_t &mat, const cu_real_t temp) {
   const cu_real_t alpha = mat.temperature_rescaling_alpha;
   const cu_real_t tc = mat.temperature_rescaling_Tc;
   const cu_real_t sigma = mat.H_th_sigma;
   #ifdef CUDA_DP
      const double resc_temp = (temp < tc) ? tc * pow(temp / tc, alpha) : temp;
      return sigma * sqrt(resc_temp);
   #else
      const float resc_temp = (temp < tc) ? tc * __powf(temp / tc, alpha) : temp;
      return sigma * sqrtf(resc_temp);
   #endif
}

// temperature of the Gaussian HAMR laser profile centred on the head (px, py) at (cx, cy)
inline __device__ cu_real_t laser_temperature(const cu_real_t cx, const cu_real_t cy, const cu_real_t px, const cu_real_t py,
                                              const cu_real_t Tmin, const cu_real_t DeltaT,
                                              const cu_real_t one_over_denx, const cu_real_t one_over_deny) {
   const cu_real_t cx2 = (cx - px) * (cx - px);
   const cu_real_t cy2 = (cy - py) * (cy - py);
   #ifdef CUDA_DP
      const cu_real_t exp_x = exp(-cx2 * one_over_denx);
      const cu_real_t exp_y = exp(-cy2 * one_over_deny);
   #else
      const cu_real_t exp_x = __expf(-cx2 * one_over_denx);
      const cu_real_t exp_y = __expf(-cy2 * one_over_deny);
   #endif
   return Tmin + DeltaT * exp_x * exp_y;
}

/*
 * Shared kernel definitions
 */

__global__ void init_rng(curandState *state, int seed);

__global__ void update_non_exchange_spin_fields_kernel(
    int *material, material_parameters_t *material_params, cu_real_t *x_spin,
    cu_real_t *y_spin, cu_real_t *z_spin, cu_real_t *x_sp_field,
    cu_real_t *y_sp_field, cu_real_t *z_sp_field, int num_atoms, int first_atom);

__global__ void update_isotropic_exchange_and_anisotropy_fields_kernel(
    const int *material, const material_parameters_t *material_params,
    const int *nbr_start, const int *nbr_atom, const cu_real_t *nbr_J,
    const cu_real_t *x_spin, const cu_real_t *y_spin, const cu_real_t *z_spin,
    cu_real_t *x_sp_field, cu_real_t *y_sp_field, cu_real_t *z_sp_field,
    int num_atoms, int first_atom, const int *active);

__global__ void
update_external_fields_kernel(cu_real_t *x_dip_field, cu_real_t *y_dip_field,
                              cu_real_t *z_dip_field, cu_real_t *x_ext_field,
                              cu_real_t *y_ext_field, cu_real_t *z_ext_field,
                              int num_atoms);

__global__ void update_cell_magnetization(
    cu_real_t *x_spin, cu_real_t *y_spin, cu_real_t *z_spin, int *material,
    int *cell, material_parameters_t *material_params, cu_real_t *x_mag,
    cu_real_t *y_mag, cu_real_t *z_mag, int num_atoms);

__global__ void update_dipolar_fields(
    cu_real_t *x_mag, cu_real_t *y_mag, cu_real_t *z_mag, cu_real_t *x_coord,
    cu_real_t *y_coord, cu_real_t *z_coord, cu_real_t *volume,
    cu_real_t *x_cell_field, cu_real_t *y_cell_field, cu_real_t *z_cell_field,
    cu_real_t *x_cell_mu0H_field, cu_real_t *y_cell_mu0H_field,
    cu_real_t *z_cell_mu0H_field, cu_real_t *d_tensor_xx,
    cu_real_t *d_tensor_xy, cu_real_t *d_tensor_xz, cu_real_t *d_tensor_yy,
    cu_real_t *d_tensor_yz, cu_real_t *d_tensor_zz, int *d_cell_id_array,
    int *d_num_atoms_in_cell, int n_local_cells, int n_cells);

__global__ void update_atomistic_dipolar_fields(
    cu_real_t *x_cell_field, cu_real_t *y_cell_field, cu_real_t *z_cell_field,
    cu_real_t *x_cell_mu0H_field, cu_real_t *y_cell_mu0H_field,
    cu_real_t *z_cell_mu0H_field, cu_real_t *x_dip_field,
    cu_real_t *y_dip_field, cu_real_t *z_dip_field, cu_real_t *x_mu0H_dip_field,
    cu_real_t *y_mu0H_dip_field, cu_real_t *z_mu0H_dip_field, int *cells,
    int n_atoms);

__global__ void apply_global_temperature_kernel(
    cu_real_t *x_field_array, cu_real_t *y_field_array,
    cu_real_t *z_field_array, cu_real_t temperature, curandState *rand_states,
    material_parameters_t *material_params, int *material, int n_atoms,
    const int *active);

__global__ void update_applied_fields_kernel(
    cu_real_t *x_field_array, cu_real_t *y_field_array,
    cu_real_t *z_field_array, const cu_real_t Hx, const cu_real_t Hy,
    const cu_real_t Hz, int *material,
    vcuda::internal::material_parameters_t *material_params, const int n_atoms);

__global__ void update_fmr_fields_kernel(
    cu_real_t *x_field_array, cu_real_t *y_field_array,
    cu_real_t *z_field_array, const cu_real_t Hx, const cu_real_t Hy,
    const cu_real_t Hz, const double real_time, const bool local_fmr_field,
    int *material, material_parameters_t *material_params, const int n_atoms);

__global__ void update_neel_fields_kernel(cu_real_t *d_neel_tensor,
                                          cu_real_t *x_spin, cu_real_t *y_spin,
                                          cu_real_t *z_spin,
                                          cu_real_t *x_sp_field,
                                          cu_real_t *y_sp_field,
                                          cu_real_t *z_sp_field, int n_atoms,
                                          int first_atom);

__global__ void update_rotational_fields_kernel(
    int *material, material_parameters_t *material_params,
    material_parameters_rotational_t *material_params_r, cu_real_t *x_spin,
    cu_real_t *y_spin, cu_real_t *z_spin, cu_real_t *x_sp_field,
    cu_real_t *y_sp_field, cu_real_t *z_sp_field, int n_atoms, int first_atom);

namespace stats {
extern bool use_cpu;
extern bool device_owns_spins;
extern bool host_spins_current;
void finalise();
} // namespace stats

} // namespace internal

#endif

} // namespace vcuda

#endif // CUDA_INTERNAL_H_
