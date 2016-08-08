//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>

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
// Host function to calculate external fields using gpu kernel
//------------------------------------------------------------------------------
void update_external_fields (){

    // Find the addresses in the device address space
   int * d_materials = thrust::raw_pointer_cast(cu::atoms::type_array.data());
   cu::material_parameters_t * d_material_params = thrust::raw_pointer_cast (cu::mp::materials.data());

   cu_real_t * d_x_dip_field = thrust::raw_pointer_cast(cu::x_dipolar_field_array.data());
   cu_real_t * d_y_dip_field = thrust::raw_pointer_cast(cu::y_dipolar_field_array.data());
   cu_real_t * d_z_dip_field = thrust::raw_pointer_cast(cu::z_dipolar_field_array.data());

   cu_real_t * d_x_ext_field = thrust::raw_pointer_cast(cu::x_total_external_field_array.data());
   cu_real_t * d_y_ext_field = thrust::raw_pointer_cast(cu::y_total_external_field_array.data());
   cu_real_t * d_z_ext_field = thrust::raw_pointer_cast(cu::z_total_external_field_array.data());

   // copy simulation variables to temporary constants
   const cu_real_t global_temperature = sim::temperature;
   const cu_real_t Hx = sim::H_vec[0]*sim::H_applied;
   const cu_real_t Hy = sim::H_vec[1]*sim::H_applied;
   const cu_real_t Hz = sim::H_vec[2]*sim::H_applied;
   const int num_atoms = ::atoms::num_atoms;

   // Call kernel to calculate external fields
   cu::update_external_fields_kernel <<< cu::grid_size, cu::block_size >>> (
         d_materials,
         d_material_params,
         d_x_dip_field, d_y_dip_field, d_z_dip_field,
         d_x_ext_field, d_y_ext_field, d_z_ext_field,
         cu::d_rand_state,
         global_temperature,
         Hx, Hy, Hz,
         num_atoms);

   // Check for errors
   check_cuda_errors (__FILE__, __LINE__);

   // update dipole field
   update_dipolar_fields();

   // std::ofstream fields("should_be_normal.txt");
   // for (size_t i = 0; i < cu::x_total_external_field_array.size(); ++i) {
   //    fields << cu::x_total_external_field_array[i] << std::endl;
   // }
   // fields.close();

   return;

}

//------------------------------------------------------------------------------
// Kernel function to calculate external fields
//------------------------------------------------------------------------------
__global__ void update_external_fields_kernel (
      int *  material,
      vcuda::internal::material_parameters_t * material_params,
      cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
      cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
      curandState * rand_states,
      cu_real_t global_temperature,
      cu_real_t Hx_app, cu_real_t Hy_app, cu_real_t Hz_app,
      int n_atoms
      )
{

   // Thread identification
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   for ( size_t atom = tid;
         atom < n_atoms;
         atom += blockDim.x * gridDim.x){

      // Get material of atom
      int mid = material[atom];
      // Load parameters to local variables from memory
      cu::material_parameters_t mat = material_params[mid];

      // initialize registers for total external field
      cu_real_t field_x = 0.0;
      cu_real_t field_y = 0.0;
      cu_real_t field_z = 0.0;

      /*
      * TODO: HAMR fields
      */

      // thermal field calculation
      //cu_real_t temp = mat.temperature;
      cu_real_t temp = global_temperature;
      cu_real_t alpha = mat.temperature_rescaling_alpha;
      cu_real_t sigma = mat.H_th_sigma;
      cu_real_t tc = mat.temperature_rescaling_Tc;

      #ifdef CUDA_DP
         double resc_temp = (temp < tc) ? tc * pow(temp / tc, alpha) : temp;
         double sq_temp = sqrt(resc_temp);
      #else
         float resc_temp = (temp < tc) ? tc * __powf(temp / tc, alpha) : temp;
         float sq_temp = sqrtf(resc_temp);
      #endif

      field_x = sigma * sq_temp * curand_normal_double (&rand_states[tid]);
      field_y = sigma * sq_temp * curand_normal_double (&rand_states[tid]);
      field_z = sigma * sq_temp * curand_normal_double (&rand_states[tid]);

      // Local applied field
      cu_real_t norm_h = mat.applied_field_strength;
      cu_real_t hx = mat.applied_field_unit_x;
      cu_real_t hy = mat.applied_field_unit_y;
      cu_real_t hz = mat.applied_field_unit_z;

      field_x += norm_h * hx;
      field_y += norm_h * hy;
      field_z += norm_h * hz;

      // Global applied field
      field_x += Hx_app;
      field_y += Hy_app;
      field_z += Hz_app;

      /*
      * TODO: FMR fields?
      */

      /*
      * Dipolar fields
      */

      field_x += x_dip_field[atom];
      field_y += y_dip_field[atom];
      field_z += z_dip_field[atom];

      // Write back to main memory
      x_ext_field[atom] = field_x;
      y_ext_field[atom] = field_y;
      z_ext_field[atom] = field_z;

   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
