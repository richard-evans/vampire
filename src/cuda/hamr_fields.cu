//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "cuda.hpp"
#include "errors.hpp"
#include "hamr.hpp"
#include "random.hpp"
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
	
		// Function to calculate temperature of atom depending on Gaussian profile
		__device__ cu_real_t calculate_gaussian_profile(
			cu_real_t * atoms_coord_x, cu_real_t * atoms_coord_y,
			int atom, cu_real_t Tmin, cu_real_t Tmax,
			cu_real_t laser_sigma_x2, cu_real_t laser_sigma_y2,
			cu_real_t px, cu_real_t py
			)
		{
			const cu_real_t cx = atoms_coord_x[atom];  
			const cu_real_t cy = atoms_coord_y[atom];  
			const cu_real_t cx2 = (cx-px)*(cx-px);
			const cu_real_t cy2 = (cy-py)*(cy-py);
			cu_real_t temperature = 0.0;

			#ifdef CUDA_DP
				double temp = Tmin + (Tmax-Tmin) * exp(-cx2/(2.0 * laser_sigma_x2)) * exp(-cy2/(2.0 * laser_sigma_y2));
				// cu_real_t sqrt_T = __dsqrt_rn( Tmin + (Tmax-Tmin) * exp(-cx2/(2.0 * laser_sigma_x2)) * exp(-cy2/(2.0 * laser_sigma_y2)) );
			#else
				float temp = Tmin + (Tmax-Tmin) * __expf(-cx2/(2.0 * laser_sigma_x2)) * __expf(-cy2/(2.0 * laser_sigma_y2));
				// cu_real_t sqrt_T = __fsqrt_rn( Tmin + (Tmax-Tmin) * __expf(-cx2/(2.0 * laser_sigma_x2)) * __expf(-cy2/(2.0 * laser_sigma_y2)) );
			#endif

			temperature = temp;

			return temperature;
		} // end calculate_gaussian_profile


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
				cu_real_t rsigma = sigma*__dsqrt_rn(rescaled_temperature);
			#else
				cu_real_t rescaled_temperature = temperature < Tc ? Tc* __powf(temperature/Tc,alpha) : temperature;
				cu_real_t rsigma = sigma*__frsqrt_rn(rescaled_temperature);
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

				// Initialise register to hold total field
				cu_real_t field_x = 0.0;
				cu_real_t field_y = 0.0;
				cu_real_t field_z = 0.0;

				// material dependent temperature rescaling
				const cu_real_t alpha = mat.temperature_rescaling_alpha;
				const cu_real_t Tc    = mat.temperature_rescaling_Tc;
				const cu_real_t sigma = mat.H_th_sigma;

				// Compute thermal field
				// __calculate_thermal_field(temperature, alpha, Tc, sigma, local_state, dev_field[0], dev_field[1], dev_field[2]);
				field_x = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);
				field_y = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);
				field_z = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);

				x_field_array[i] = field_x;
				y_field_array[i] = field_y;
				z_field_array[i] = field_z;
			}
			return;
		} // end apply_global_temperature_kernel


		// Kernel to apply thermal field with local temperature
		__global__ void apply_local_temperature_kernel(
			 	cu_real_t * atoms_coord_x, cu_real_t * atoms_coord_y,
				cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
				cu_real_t Tmin, cu_real_t Tmax,
				cu_real_t global_temperature,
				cu_real_t laser_sigma_x2, cu_real_t laser_sigma_y2,
				cu_real_t px, cu_real_t py,	 
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

				// Initialise register to hold total field
				cu_real_t field_x = 0.0;
				cu_real_t field_y = 0.0;
				cu_real_t field_z = 0.0;

				cu_real_t temperature = global_temperature;

				// material dependent temperature rescaling
				const cu_real_t alpha = mat.temperature_rescaling_alpha;
				const cu_real_t Tc    = mat.temperature_rescaling_Tc;
				const cu_real_t sigma = mat.H_th_sigma;

				// Assign tempeerature to atoms according to Gaussian profile
				temperature = calculate_gaussian_profile(atoms_coord_x, atoms_coord_y, i, Tmin, Tmax, laser_sigma_x2, laser_sigma_y2, px, py);

				// Compute thermal field
				// __calculate_thermal_field(temperature, alpha, Tc, sigma, local_state, field_x, field_y, field_z);
				field_x = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);
				field_y = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);
				field_z = calculate_thermal_field(temperature, alpha, Tc, sigma, local_state);

				// No need to sum over previous values because already zeroed initially
				x_field_array[i] = field_x;
				y_field_array[i] = field_y;
				z_field_array[i] = field_z;
			}

			return;
		} // end apply_local_temperature_kernel



		// Kernel to apply thermal field with local temperature
		__global__ void apply_local_external_field_kernel(
			 	cu_real_t * atoms_coord_x, cu_real_t * atoms_coord_y,
				cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
				cu_real_t Hx_app, cu_real_t Hy_app, cu_real_t Hz_app, 
				cu_real_t head_position_x, cu_real_t head_position_y,
				cu_real_t H_bounds_x, cu_real_t H_bounds_y, 
				cu_real_t NPS,
				int n_atoms
				)
		{
			int tid = blockIdx.x * blockDim.x + threadIdx.x;
			for (int i = tid; 
				i < n_atoms; 
				i += blockDim.x * gridDim.x){

				// Initialise register to hold total field
				cu_real_t field_x = 0.0;
				cu_real_t field_y = 0.0;
				cu_real_t field_z = 0.0;

				const cu_real_t cx = atoms_coord_x[i]; 
				const cu_real_t cy = atoms_coord_y[i]; 
				const cu_real_t Hloc_min_x = head_position_x - H_bounds_x - NPS;  // Shift field box in downtrack of NPS
				const cu_real_t Hloc_max_x = head_position_x + H_bounds_x - NPS;  // Shift field box in downtrack of NPS
				const cu_real_t Hloc_min_y = head_position_y - H_bounds_y;
				const cu_real_t Hloc_max_y = head_position_y + H_bounds_y;

				// If atoms within field box, add contribution from external field
				if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
					field_x = Hx_app;
					field_y = Hy_app;
					field_z = Hz_app;
				}

				x_field_array[i] += field_x;
				y_field_array[i] += field_y;
				z_field_array[i] += field_z;
			}

			return;
		} // end apply_local_external_field_kernel


		// Function to calculate update hamr field
		void update_hamr_field(
			cu_real_t temperature,
			cu_real_t Tmin, cu_real_t Tmax,
			cu_real_t Hx_app, cu_real_t Hy_app, cu_real_t Hz_app,
			int num_atoms
			)
		{
			// Check that hamr field calculation has been called
			if(err::check==true){ std::cout << "calculate_hamr_fields has been called" << std::endl;}

			// Initialise to zero hamr fields
			cudaMemset(cu::d_x_hamr_field, 0, num_atoms * sizeof(cu_real_t));
			cudaMemset(cu::d_y_hamr_field, 0, num_atoms * sizeof(cu_real_t));
			cudaMemset(cu::d_z_hamr_field, 0, num_atoms * sizeof(cu_real_t));

			check_cuda_errors (__FILE__, __LINE__);

			if(::hamr::head_laser_on){

				check_cuda_errors (__FILE__, __LINE__);

				// Determine constants
				const cu_real_t laser_sigma_x2 = cu::hamr::d_laser_sigma_x * cu::hamr::d_laser_sigma_x;
				const cu_real_t laser_sigma_y2 = cu::hamr::d_laser_sigma_y * cu::hamr::d_laser_sigma_y;
				const cu_real_t px = cu::hamr::d_head_position_x;
				const cu_real_t py = cu::hamr::d_head_position_y;

				// Apply thermal field
				apply_local_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::atoms::d_x_coord, cu::atoms::d_x_coord,
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					Tmin, Tmax,
					temperature,
					laser_sigma_x2, laser_sigma_y2,
					px, py,	 
					cu::d_rand_state,
					cu::mp::d_material_params,
					cu::atoms::d_materials,
					num_atoms);

				check_cuda_errors (__FILE__, __LINE__);

				// Apply external field 
				apply_local_external_field_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::atoms::d_x_coord, cu::atoms::d_x_coord,
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					Hx_app, Hy_app, Hz_app,
					cu::hamr::d_head_position_x, cu::hamr::d_head_position_y,
					cu::hamr::d_H_bounds_x, cu::hamr::d_H_bounds_y, 
					cu::hamr::d_NPS,
					num_atoms);

				check_cuda_errors (__FILE__, __LINE__);
			}
			// Apply global temperature if laser is off
			else{

				check_cuda_errors (__FILE__, __LINE__);

				apply_global_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					temperature,
					cu::d_rand_state,
					cu::mp::d_material_params,
					cu::atoms::d_materials,
					num_atoms);

				check_cuda_errors (__FILE__, __LINE__);
			}

			return;
		} // end update_hamr_field


	} // end namespace internal

} // end namespace vcuda


#endif