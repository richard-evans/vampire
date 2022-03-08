//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
// Reviewed: Andrea Meo 2022
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
// Host function to calculate spin fields using gpu kernel
//------------------------------------------------------------------------------
void update_spin_fields ()
{

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
// Device function for uniaxial anisotropy
//------------------------------------------------------------------------------

__device__ cu_real_t uniaxial_anisotropy_energy(cu::material_parameters_t &material,
        cu_real_t sx, cu_real_t sy, cu_real_t sz)
{
      // Factors for k4
      const cu_real_t fiveothirtyfive  = 0.14285714285; // 5/35
      const cu_real_t thirtyothirtyfive = 0.85714285714; // 30/35

      const cu_real_t ex = material.anisotropy_unit_x;
      const cu_real_t ey = material.anisotropy_unit_y;
      const cu_real_t ez = material.anisotropy_unit_z;

      const cu_real_t sdote = sx * ex + sy * ey + sz * ez;
      const cu_real_t sdote2 = sdote * sdote;

      const cu_real_t k2 = material.sh2;  
      const cu_real_t k4 = material.sh4;
      const cu_real_t k6 = material.sh6;

      const cu_real_t Ek2 = -k2 * sdote2;
      const cu_real_t Ek4 =  k4 * (sdote2*sdote2 - thirtyothirtyfive*sdote2 - fiveothirtyfive);
      const cu_real_t Ek6 = -0.04166666666 * k6 * (231.0*sdote2*sdote2*sdote2 - 315.0*sdote2*sdote2 + 105.0*sdote2); // factor = 2/3 * -1/16 = -1/6 = -0.04166666666

      cu_real_t E = material.mu_s_si * (Ek2 + Ek4 + Ek6);

      return E;
}





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

      const cu_real_t sx = x_spin[i];
      const cu_real_t sy = y_spin[i];
      const cu_real_t sz = z_spin[i];

      const cu_real_t ex = material.anisotropy_unit_x;
      const cu_real_t ey = material.anisotropy_unit_y;
      const cu_real_t ez = material.anisotropy_unit_z;

      const cu_real_t sdote = sx * ex + sy * ey + sz * ez;
      const cu_real_t sdote3 = sdote * sdote * sdote;
      const cu_real_t sdote5 = sdote3 * sdote * sdote;

      /*
       * Spherical harmonics
       */

      const cu_real_t scale = 0.6666666666666667;

	// Reduced anisotropy constants ku/mu_s [J/T]
      const cu_real_t k2 = material.sh2;
      const cu_real_t k4 = material.sh4;
      const cu_real_t k6 = material.sh6;

      const cu_real_t ek2 = k2 * 3.0 * sdote;
      const cu_real_t ek4 = -k4 * 0.125 * (140.0 * sdote3 - 60.0 *sdote);
      const cu_real_t ek6 = k6 * 0.0625 * (1386.0 * sdote5 - 1260.0 * sdote3 + 210.0 * sdote);

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
