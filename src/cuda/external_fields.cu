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
#include "thermal_fields.hpp"

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
   //int * d_materials = thrust::raw_pointer_cast(cu::atoms::type_array.data());
   //cu::material_parameters_t * d_material_params = thrust::raw_pointer_cast (cu::mp::materials.data());

   /*
   cu_real_t * d_x_dip_field = thrust::raw_pointer_cast(cu::x_dipolar_field_array.data());
   cu_real_t * d_y_dip_field = thrust::raw_pointer_cast(cu::y_dipolar_field_array.data());
   cu_real_t * d_z_dip_field = thrust::raw_pointer_cast(cu::z_dipolar_field_array.data());

   cu_real_t * d_x_ext_field = thrust::raw_pointer_cast(cu::x_total_external_field_array.data());
   cu_real_t * d_y_ext_field = thrust::raw_pointer_cast(cu::y_total_external_field_array.data());
   cu_real_t * d_z_ext_field = thrust::raw_pointer_cast(cu::z_total_external_field_array.data());
   */
   // copy simulation variables to temporary constants
   const cu_real_t Hx = sim::H_vec[0]*sim::H_applied;
   const cu_real_t Hy = sim::H_vec[1]*sim::H_applied;
   const cu_real_t Hz = sim::H_vec[2]*sim::H_applied;
   const int num_atoms = ::atoms::num_atoms;

   // update hamr fields (Happ and Hth) if program is hamr simulations
	if(sim::program==7){
      cu::update_hamr_field();
   }
   // Calculate global thermal field otherwise
   else{
      cu::update_global_thermal_field();
   }

//   // update dipole field
//    update_dipolar_fields();  //-- disabled  as causes NaN and deferred to CPU code for now

   // Call kernel to calculate external fields
   cu::update_external_fields_kernel <<< cu::grid_size, cu::block_size >>> (
         cu::atoms::d_materials, cu::mp::d_material_params,
         cu::d_x_dip_field, cu::d_y_dip_field, cu::d_z_dip_field,
         cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
         cu::d_x_thermal_field, cu::d_y_thermal_field, cu::d_z_thermal_field,
         cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
         Hx, Hy, Hz,
         num_atoms);

   // Check for errors
   check_cuda_errors (__FILE__, __LINE__);

   // std::ofstream fields("should_be_normal.txt");
   // for (size_t i = 0; i < cu::x_total_external_field_array.size(); ++i) {
   //    fields << cu::x_total_external_field_array[i] << std::endl;
   // }
   // fields.close();

   //std::cerr << num_atoms << "\t";
   //std::cerr << cu::x_total_external_field_array[0] << "\t";
   //std::cerr << cu::y_total_external_field_array[0] << "\t";
   //std::cerr << cu::z_total_external_field_array[0] << "\t";


   return;

}

//------------------------------------------------------------------------------
// Kernel function to calculate external fields
//------------------------------------------------------------------------------
__global__ void update_external_fields_kernel (
      int *  material,
      vcuda::internal::material_parameters_t * material_params,
      cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
      cu_real_t * x_hamr_field, cu_real_t * y_hamr_field, cu_real_t * z_hamr_field,
      cu_real_t * x_thermal_field, cu_real_t * y_thermal_field, cu_real_t * z_thermal_field,
      cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
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

      // Hamr fields
      field_x += x_hamr_field[atom];
      field_y += y_hamr_field[atom];
      field_z += z_hamr_field[atom];
//      printf("       %d  %lf  %lf  %lf\n",atom,x_hamr_field[atom],y_hamr_field[atom],z_hamr_field[atom]);

      // Thermal fields
      field_x += x_thermal_field[atom];
      field_y += y_thermal_field[atom];
      field_z += z_thermal_field[atom];
      // if(atom<10){printf("       %d  %lf  %lf  %lf\n",atom,x_thermal_field[atom],y_thermal_field[atom],z_thermal_field[atom]);}

      // Local applied field
      cu_real_t norm_h = mat.applied_field_strength;

      field_x += norm_h * mat.applied_field_unit_x;
      field_y += norm_h * mat.applied_field_unit_y;
      field_z += norm_h * mat.applied_field_unit_z;

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
      
//      printf("       %d  %lf  %lf  %lf\n",atom,x_dip_field[atom],y_dip_field[atom],z_dip_field[atom]);

      // Write back to main memory
      x_ext_field[atom] = field_x;
      y_ext_field[atom] = field_y;
      z_ext_field[atom] = field_z;

   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
