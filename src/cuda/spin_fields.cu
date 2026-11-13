//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015.
//   All rights reserved.
//
//   Reviewed: Andrea Meo 2022
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "exchange_fields.hpp"
#include "data.hpp"
#include "internal.hpp"

#include "spin_fields.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

//------------------------------------------------------------------------------
// Host function to calculate spin fields of atoms first_atom to num_atoms-1
// using gpu kernels (the exchange matrix path always computes all atoms)
//------------------------------------------------------------------------------
void update_spin_fields (const int num_atoms, const int first_atom)
{

   if(!cu::exchange::exchange_initialised) cu::exchange::initialise_exchange();

   // Isotropic exchange: exchange and anisotropy fields in one kernel, which
   // overwrites the spin field
   if(cu::exchange::use_isotropic_kernel){
      cu::update_isotropic_exchange_and_anisotropy_fields_kernel <<< cu::grid_size, cu::block_size >>> (
            cu::atoms::d_materials, cu::mp::d_material_params,
            cu::exchange::d_nbr_start, cu::exchange::d_nbr_atom, cu::exchange::d_nbr_J,
            cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
            cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
            cu::num_integrated_atoms(num_atoms), cu::num_integrated_atoms(first_atom), cu::d_active_atoms);

      check_cuda_errors (__FILE__, __LINE__);

      //If neel fields are enabled, execute the kernel
      if (anisotropy::is_neel_enabled()) cu::update_neel_fields(num_atoms, first_atom);

      //If rotational fields are enabled, execute the kernel
      if (anisotropy::is_rotational_enabled()) cu::update_rotational_fields(num_atoms, first_atom);

      return;
   }

   // This exchange field zero out stuff, so they should come first
   cu::exchange::calculate_exchange_fields ();

   check_cuda_errors (__FILE__, __LINE__);

   //If neel fields are enabled, execute the kernel
   if (anisotropy::is_neel_enabled()) cu::update_neel_fields(num_atoms, first_atom);

   //If rotational fields are enabled, execute the kernel
   if (anisotropy::is_rotational_enabled()) cu::update_rotational_fields(num_atoms, first_atom);

   // Call kernel to calculate non-exchange spin fields
   cu::update_non_exchange_spin_fields_kernel <<< cu::grid_size, cu::block_size >>> (
         cu::atoms::d_materials, cu::mp::d_material_params,
         cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
         cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
         num_atoms, first_atom);

   check_cuda_errors (__FILE__, __LINE__);

}

//------------------------------------------------------------------------------
// Device function for uniaxial anisotropy
//------------------------------------------------------------------------------

__device__ cu_real_t uniaxial_anisotropy_energy(cu::material_parameters_t &material,
        cu_real_t sx, cu_real_t sy, cu_real_t sz)
{
      // same expressions as the CPU (anisotropy/uniaxial_order_2, _4, _6)
      const cu_real_t six_o_seven = 6.0 / 7.0;
      const cu_real_t fifteen_o_eleven = 15.0 / 11.0;
      const cu_real_t five_o_eleven = 5.0 / 11.0;

      const cu_real_t ex = material.anisotropy_unit_x;
      const cu_real_t ey = material.anisotropy_unit_y;
      const cu_real_t ez = material.anisotropy_unit_z;

      const cu_real_t sdote = sx * ex + sy * ey + sz * ez;
      const cu_real_t sdote2 = sdote * sdote;

      const cu_real_t k2 = material.sh2;
      const cu_real_t k4 = material.sh4;
      const cu_real_t k6 = material.sh6;

      const cu_real_t sdote4 = sdote2 * sdote2;
      const cu_real_t Ek2 = -k2 * sdote2;
      const cu_real_t Ek4 = -k4 * (sdote4 - six_o_seven*sdote2);
      const cu_real_t Ek6 = -k6 * (sdote2*sdote4 - fifteen_o_eleven*sdote4 + five_o_eleven*sdote2);

      cu_real_t E = material.mu_s_si * (Ek2 + Ek4 + Ek6);

      return E;
}

//------------------------------------------------------------------------------
// Device function to calculate the uniaxial and cubic anisotropy field of a
// spin
//------------------------------------------------------------------------------
__device__ __forceinline__ void anisotropy_field (
      const cu::material_parameters_t & material,
      const cu_real_t sx, const cu_real_t sy, const cu_real_t sz,
      cu_real_t & field_x, cu_real_t & field_y, cu_real_t & field_z
      )
{

      field_x = 0.0;
      field_y = 0.0;
      field_z = 0.0;

      const cu_real_t ex = material.anisotropy_unit_x;
      const cu_real_t ey = material.anisotropy_unit_y;
      const cu_real_t ez = material.anisotropy_unit_z;

      const cu_real_t sdote = sx * ex + sy * ey + sz * ez;
      const cu_real_t sdote3 = sdote * sdote * sdote;
      const cu_real_t sdote5 = sdote3 * sdote * sdote;

      /*
       * Spherical harmonics
       */

      const cu_real_t scale = 2.0 / 3.0;

   // Reduced anisotropy constants ku/mu_s [J/T]
      const cu_real_t k2 = material.sh2;
      const cu_real_t k4 = material.sh4;
      const cu_real_t k6 = material.sh6;

      const cu_real_t ek2 = k2 * 3.0 * sdote;
      const cu_real_t ek4 = k4 * (6.0 * sdote3 - ((18.0/7.0) * sdote));
      const cu_real_t ek6 = k6 * (3.0/11.0) * (((33.0)*sdote5) - (30.0 * sdote3) + (5.0*sdote));

      field_x += scale * ex * (ek2 + ek4 + ek6);
      field_y += scale * ey * (ek2 + ek4 + ek6);
      field_z += scale * ez * (ek2 + ek4 + ek6);

      /*
       * Cubic Anisotropy
       */

      // scale factor from derivative of E = -1/2 (sx^4 + sy^4 + sz^4)
      const cu_real_t scalec = 0.5*4.0;

      const cu_real_t kc4 = material.kc4;
      const cu_real_t kc6 = material.kc6;

      const cu_real_t sx2 = sx*sx;
      const cu_real_t sy2 = sy*sy;
      const cu_real_t sz2 = sz*sz;

      field_x += scalec * sx * ((sx2 * kc4) + (-sy2 * sz2 * kc6));
      field_y += scalec * sy * ((sy2 * kc4) + (-sz2 * sx2 * kc6));
      field_z += scalec * sz * ((sz2 * kc4) + (-sx2 * sy2 * kc6));

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

}

__global__ void update_non_exchange_spin_fields_kernel (
      int * material,
      cu::material_parameters_t * material_params,
      cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
      cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
      int n_atoms, int first_atom
      )
{
   for ( int i = first_atom + blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x){

      cu_real_t field_x, field_y, field_z;
      anisotropy_field(material_params[material[i]], x_spin[i], y_spin[i], z_spin[i], field_x, field_y, field_z);

      x_sp_field[i] += field_x;
      y_sp_field[i] += field_y;
      z_sp_field[i] += field_z;

   }
}

//------------------------------------------------------------------------------
// Kernel calculating isotropic exchange and anisotropy fields in one pass,
// overwriting the spin field. Neighbours of each atom are summed in increasing
// atom order, as in the rows of the exchange matrix.
//------------------------------------------------------------------------------
__global__ void update_isotropic_exchange_and_anisotropy_fields_kernel (
      const int * material, const cu::material_parameters_t * material_params,
      const int * nbr_start, const int * nbr_atom, const cu_real_t * nbr_J,
      const cu_real_t * x_spin, const cu_real_t * y_spin, const cu_real_t * z_spin,
      cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
      int n_atoms, int first_atom, const int * active
      )
{
   for ( int k = first_atom + blockIdx.x * blockDim.x + threadIdx.x;
         k < n_atoms;
         k += blockDim.x * gridDim.x){

      const int i = active ? active[k] : k;

      // atoms that are not integrated need no field
      if(material_params[material[i]].skip) continue;

      // exchange field
      cu_real_t hx = 0.0;
      cu_real_t hy = 0.0;
      cu_real_t hz = 0.0;
      const int end = nbr_start[i + 1];
      for(int n = nbr_start[i]; n < end; n++){
         const int j = nbr_atom[n];
         const cu_real_t J = nbr_J[n];
         hx += J * x_spin[j];
         hy += J * y_spin[j];
         hz += J * z_spin[j];
      }

      cu_real_t field_x, field_y, field_z;
      anisotropy_field(material_params[material[i]], x_spin[i], y_spin[i], z_spin[i], field_x, field_y, field_z);

      x_sp_field[i] = hx + field_x;
      y_sp_field[i] = hy + field_y;
      z_sp_field[i] = hz + field_z;

   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
