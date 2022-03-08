//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>

// Vampire headers
#include "cuda.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "typedefs.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA
// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

	// module internal namespace
	namespace internal{
	
		// Function to calculate thermal field, also accounting for rescaling, once temperature is given
		__device__ cu_real_t calculate_thermal_field(
			const cu_real_t temperature, 
			const cu_real_t alpha, const cu_real_t Tc, const cu_real_t sigma,
			curandState local_state
			)
		{
			cu_real_t field = 0.0;

			// Determine sigma accounting for rescaling if any
			#ifdef CUDA_DP
				cu_real_t rescaled_temperature = temperature < Tc ? Tc*pow(temperature/Tc,alpha) : temperature;
				cu_real_t rsigma = sigma*sqrt(rescaled_temperature);
			#else
				cu_real_t rescaled_temperature = temperature < Tc ? Tc* __powf(temperature/Tc,alpha) : temperature;
				cu_real_t rsigma = sigma*sqrtf(rescaled_temperature);
			#endif

			// Determine field 
			#ifdef CUDA_DP
				double field_tmp = rsigma * curand_normal_double (&local_state);
			#else
			   float field_tmp = rsigma * curand_normal(&local_state);
			#endif

			field = field_tmp;
			return field;
		} // end  __calculate_thermal_field


		// Kernel to apply thermal field with global temperature
		__global__ void apply_global_temperature_kernel(
				cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
				cu_real_t temperature,
				curandState * rand_states,
				cu::material_parameters_t * material_params,
				int * material,
				int n_atoms
				)
		{
			int tid = blockIdx.x * blockDim.x + threadIdx.x;
			for (int i = tid; 
				i < n_atoms; 
				i += blockDim.x * gridDim.x){

				// Get material of atom i
				int mid = material[i];
				// Load parameters from memory
				cu::material_parameters_t mat = material_params[mid];

				// Load the curand state into local memory
				curandState local_state = rand_states[tid];

				cu_real_t field_x = 0.0;
				cu_real_t field_y = 0.0;
				cu_real_t field_z = 0.0;

      		// thermal field calculation
      		cu_real_t temp = temperature;
      		cu_real_t alpha = mat.temperature_rescaling_alpha;
      		cu_real_t sigma = mat.H_th_sigma;
      		cu_real_t tc = mat.temperature_rescaling_Tc;

      		#ifdef CUDA_DP
      		   double resc_temp = (temp < tc) ? tc * pow(temp / tc, alpha) : temp;
      		   double rsigma = sigma*sqrt(resc_temp);
      		#else
      		   float resc_temp = (temp < tc) ? tc * __powf(temp / tc, alpha) : temp;
      		   float rsigma = sigma*sqrtf(resc_temp);
      		#endif

      		#ifdef CUDA_DP
      		   field_x = rsigma * curand_normal_double (&local_state);
      		   field_y = rsigma * curand_normal_double (&local_state);
      		   field_z = rsigma * curand_normal_double (&local_state);
      		#else
      		   field_x = rsigma * curand_normal(&local_state);
      		   field_y = rsigma * curand_normal(&local_state);
      		   field_z = rsigma * curand_normal(&local_state);
      		#endif
      		// if(i<10){printf("         %d  %lf  %lf  %lf\n",i,field_x,field_y,field_z);}

      		x_field_array[i] = field_x;
      		y_field_array[i] = field_y;
      		z_field_array[i] = field_z;

      		// Write local curand state back to global memory
      		rand_states[tid] = local_state;

			}
			return;
		} // end apply_global_temperature_kernel


		// Function to calculate global thermal field
		void update_global_thermal_field()
		{
			// Check that hamr field calculation has been called
			if(err::check==true){ std::cout << "update_global_thermal_field has been called" << std::endl;}

			// copy simulation variables to temporary constants
			const cu_real_t global_temperature = sim::temperature;
			const int num_atoms = ::atoms::num_atoms;

			apply_global_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
				cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
				global_temperature,
				cu::d_rand_state,
				cu::mp::d_material_params,
				cu::atoms::d_materials,
				num_atoms); 

			check_cuda_errors (__FILE__, __LINE__);

			return;
		} // end update_global_thermal_field


	} // end namespace internal

} // end namespace vcuda


#endif