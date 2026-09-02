//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
// Reviewd: Andrea Meo 2022
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
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
   // Kernel function to calculate external applied fields
   //------------------------------------------------------------------------------
   __global__ void update_applied_fields_kernel(
   		cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
         const cu_real_t Hx_app, const cu_real_t Hy_app, const cu_real_t Hz_app,
         int *  material, vcuda::internal::material_parameters_t * material_params,
   		const int n_atoms
   		)
   {
   	for ( int atom = blockIdx.x * blockDim.x + threadIdx.x;
   	      atom < n_atoms;
   	      atom += blockDim.x * gridDim.x)
   	{
         // Get material of atom
         int mid = material[atom];
         // Load parameters to local variables from memory
         cu::material_parameters_t mat = material_params[mid];
 
         // initialize registers for total external field
         cu_real_t field_x = 0.0;
         cu_real_t field_y = 0.0;
         cu_real_t field_z = 0.0;

         // Local applied field
         cu_real_t norm_h = mat.applied_field_strength;
         field_x += norm_h * mat.applied_field_unit_x;
         field_y += norm_h * mat.applied_field_unit_y;
         field_z += norm_h * mat.applied_field_unit_z;

         // Global applied field
         field_x += Hx_app;
         field_y += Hy_app;
         field_z += Hz_app;

   	   x_field_array[atom] += field_x;
   	   y_field_array[atom] += field_y;
   	   z_field_array[atom] += field_z;
   	}
   } // end calculate_applied_fields_kernel


   // Host function to calculate applied fields
   void update_applied_fields()
   {

   	// copy simulation variables to temporary constants
      const cu_real_t Hx_app = sim::H_vec[0]*sim::H_applied;
      const cu_real_t Hy_app = sim::H_vec[1]*sim::H_applied;
      const cu_real_t Hz_app = sim::H_vec[2]*sim::H_applied;

   	const int num_atoms = ::atoms::num_atoms;

   	update_applied_fields_kernel <<< cu::grid_size, cu::block_size >>> (
   		cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
         Hx_app, Hy_app, Hz_app,
         cu::atoms::d_materials, cu::mp::d_material_params,
   		num_atoms); 

   	check_cuda_errors (__FILE__, __LINE__);

   	return;
   } // end update_global_thermal_field


   //------------------------------------------------------------------------------
   // Kernel function to add the global and material-specific FMR fields
   //------------------------------------------------------------------------------
   __global__ void update_fmr_fields_kernel(
         cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
         const cu_real_t Hx_fmr, const cu_real_t Hy_fmr, const cu_real_t Hz_fmr,
         const double real_time, const bool local_fmr_field,
         int * material, material_parameters_t * material_params,
         const int n_atoms)
   {
      for (int atom = blockIdx.x * blockDim.x + threadIdx.x;
           atom < n_atoms;
           atom += blockDim.x * gridDim.x)
      {
         cu_real_t field_x = Hx_fmr;
         cu_real_t field_y = Hy_fmr;
         cu_real_t field_z = Hz_fmr;

         if (local_fmr_field)
         {
            const int mid = material[atom];
            const material_parameters_t mat = material_params[mid];
            const cu_real_t phase = cu_real_t(2.0 * M_PI) * mat.fmr_field_frequency * real_time;

            #ifdef CUDA_DP
               const cu_real_t Hsinwt = mat.fmr_field_strength * sin(phase);
            #else
               const cu_real_t Hsinwt = mat.fmr_field_strength * sinf(phase);
            #endif

            field_x += Hsinwt * mat.fmr_field_unit_x;
            field_y += Hsinwt * mat.fmr_field_unit_y;
            field_z += Hsinwt * mat.fmr_field_unit_z;
         }

         x_field_array[atom] += field_x;
         y_field_array[atom] += field_y;
         z_field_array[atom] += field_z;
      }
   }


   // Host function to calculate and add the FMR field
   void update_fmr_fields()
   {
      const double real_time = sim::time * ::mp::dt_SI;
      const double phase = (2.0 * M_PI) *
                              sim::fmr_field_frequency * real_time;
      const double Hsinwt = sim::fmr_field_strength * sin(phase);

      // Update CPU up to date for output.
      sim::fmr_field = Hsinwt;

      const cu_real_t Hx_fmr = sim::fmr_field_unit_vector[0] * Hsinwt;
      const cu_real_t Hy_fmr = sim::fmr_field_unit_vector[1] * Hsinwt;
      const cu_real_t Hz_fmr = sim::fmr_field_unit_vector[2] * Hsinwt;

      update_fmr_fields_kernel <<< cu::grid_size, cu::block_size >>> (
            cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
            Hx_fmr, Hy_fmr, Hz_fmr,
            real_time, sim::local_fmr_field,
            cu::atoms::d_materials, cu::mp::d_material_params,
            ::atoms::num_atoms);

      check_cuda_errors (__FILE__, __LINE__);
   }


   //------------------------------------------------------------------------------
   // Kernel function to calculate external fields
   //------------------------------------------------------------------------------
   __global__ void update_external_fields_kernel (
         cu_real_t * x_dip_field, cu_real_t * y_dip_field, cu_real_t * z_dip_field,
         cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
         int n_atoms
         )
   {

      // Thread identification
      int tid = blockIdx.x * blockDim.x + threadIdx.x;
      for ( size_t atom = tid;
            atom < n_atoms;
            atom += blockDim.x * gridDim.x){

         /*
         * TODO: FMR fields?
         */

			// Add dipolar fields
         x_ext_field[atom] += x_dip_field[atom];
         y_ext_field[atom] += y_dip_field[atom];
         z_ext_field[atom] += z_dip_field[atom];

      }
   }

   //------------------------------------------------------------------------------
   // Host function to calculate external fields using gpu kernel
   //------------------------------------------------------------------------------
   void update_external_fields (){

      // copy simulation variables to temporary constants
      const int num_atoms = ::atoms::num_atoms;

      // update hamr fields (Happ and Hth) if program is hamr simulations
      if(program::program==7){
         cu::update_hamr_field();
      }
      // Otherwise calculate global thermal field and applied fields
      else{
         cu::update_global_thermal_field();
         cu::update_applied_fields();
      }

      // Add the oscillating field only after it has been enabled by the FMR program.
      if(sim::enable_fmr) cu::update_fmr_fields();

   //   // update dipole field
   //    update_dipolar_fields();  //-- disabled  as causes NaN and deferred to CPU code for now

      // Call kernel to calculate external fields
      cu::update_external_fields_kernel <<< cu::grid_size, cu::block_size >>> (
            cu::d_x_dip_field, cu::d_y_dip_field, cu::d_z_dip_field,
            cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
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


} // end of internal namespace

} // end of vcuda namespace

#endif
