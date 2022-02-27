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
#include "sim.hpp"
#include "vio.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "thermal_fields.hpp"
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
			const int atom, const cu_real_t cx, const cu_real_t cy,
			const cu_real_t px, const cu_real_t py, 
			const cu_real_t Tmin, const cu_real_t DeltaT,
			const cu_real_t laser_sigma_x2, const cu_real_t laser_sigma_y2
			)
		{

			const cu_real_t cx2 = (cx-px)*(cx-px);
			const cu_real_t cy2 = (cy-py)*(cy-py);
			const cu_real_t denx = 2.0 * laser_sigma_x2;
			const cu_real_t one_over_denx = 1.0/denx;
			const cu_real_t deny = 2.0 * laser_sigma_y2;
			const cu_real_t one_over_deny = 1.0/deny;

			#ifdef CUDA_DP
				cu_real_t exp_x =  exp(-cx2*one_over_denx); 
				cu_real_t exp_y =  exp(-cy2*one_over_deny); 
			#else
				cu_real_t exp_x =  __expf(-cx2*one_over_denx); 
				cu_real_t exp_y =  __expf(-cy2*one_over_deny); 
			#endif
			// if(i<5){printf("  %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			// 						i, cx, cy, px, py, DeltaT, denx, deny, one_over_denx, one_over_deny, exp_x, exp_y); }

			cu_real_t temperature = Tmin + DeltaT * exp_x * exp_y;

			return temperature;
		} // end calculate_gaussian_profile


		// Kernel to apply thermal field with local temperature
		__global__ void apply_local_temperature_kernel(
			 	cu_real_t * atoms_coord_x, cu_real_t * atoms_coord_y,
				cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
				const cu_real_t Tmin, const cu_real_t Tmax,
				const cu_real_t laser_sigma_x2, const cu_real_t laser_sigma_y2,
				const cu_real_t px, const cu_real_t py,	 
				curandState * rand_states,
				cu::material_parameters_t * material_params,
				int * material,
				const int n_atoms
				)
		{
			int tid = blockIdx.x * blockDim.x + threadIdx.x;
			for (int i = tid; 
				i < n_atoms; 
				i += blockDim.x * gridDim.x)
			{

				// Get material of atom i
				int mid = material[i];
				// Load parameters from memory
				cu::material_parameters_t mat = material_params[mid];
				// Load the curand state into local memory
				curandState local_state = rand_states[tid];

				// Define temporary variables for field
				cu_real_t field_x = 0.0;
				cu_real_t field_y = 0.0;
				cu_real_t field_z = 0.0;

				// Assign tempeerature to atoms according to Gaussian profile
				const cu_real_t DeltaT = Tmax - Tmin;
				cu_real_t temp = calculate_gaussian_profile(i, atoms_coord_x[i], atoms_coord_y[i], px, py,
																Tmin, DeltaT, laser_sigma_x2, laser_sigma_y2);
				
				// material dependent temperature rescaling
				const cu_real_t alpha = mat.temperature_rescaling_alpha;
				const cu_real_t Tc    = mat.temperature_rescaling_Tc;
				const cu_real_t sigma = mat.H_th_sigma;

      		// thermal field calculation

      		#ifdef CUDA_DP
      		   double resc_temp = (temp < Tc) ? Tc * pow(temp / Tc, alpha) : temp;
      		   double rsigma = sigma*sqrt(resc_temp);
      		#else
      		   float resc_temp = (temp < Tc) ? Tc * __powf(temp / Tc, alpha) : temp;
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

      		x_field_array[i] = field_x;
      		y_field_array[i] = field_y;
      		z_field_array[i] = field_z;

      		// Write local curand state back to global memory
      		rand_states[tid] = local_state;
			}

			return;
		} // end apply_local_temperature_kernel



		// Kernel to apply thermal field with local temperature
		__global__ void apply_local_external_field_kernel(
			 	cu_real_t * atoms_coord_x, cu_real_t * atoms_coord_y,
				cu_real_t * x_field_array, cu_real_t * y_field_array, cu_real_t * z_field_array,
				const cu_real_t Hx_app, const cu_real_t Hy_app, const cu_real_t Hz_app, 
				const cu_real_t head_position_x, const cu_real_t head_position_y,
				const cu_real_t H_bounds_x, const cu_real_t H_bounds_y, 
				const cu_real_t NPS,
				const int n_atoms
				)
		{
			for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
				i < n_atoms; 
				i += blockDim.x * gridDim.x)
			{

				const cu_real_t cx = atoms_coord_x[i]; 
				const cu_real_t cy = atoms_coord_y[i]; 
				const cu_real_t Hloc_min_x = head_position_x - H_bounds_x - NPS;
				const cu_real_t Hloc_max_x = head_position_x + H_bounds_x - NPS;
				const cu_real_t Hloc_min_y = head_position_y - H_bounds_y;
				const cu_real_t Hloc_max_y = head_position_y + H_bounds_y;

				// If atoms within field box, add contribution from external field
				if((cx >= Hloc_min_x) && (cx <= Hloc_max_x) && (cy >= Hloc_min_y) && (cy <= Hloc_max_y)){
					x_field_array[i] += Hx_app;
					y_field_array[i] += Hy_app;
					z_field_array[i] += Hz_app;
				}

			}

			return;
		} // end apply_local_external_field_kernel


		// Function to calculate update hamr field
		void update_hamr_field(
			)
		{
			// Check that hamr field calculation has been called
			if(err::check==true){ std::cout << "calculate_hamr_fields has been called" << std::endl;}

			// copy simulation variables to temporary constants
			const cu_real_t global_temperature = sim::temperature;
			const int num_atoms = ::atoms::num_atoms;

			// // Initialise to zero hamr fields
			// cudaMemset(cu::d_x_hamr_field, 0, num_atoms * sizeof(cu_real_t));
			// cudaMemset(cu::d_y_hamr_field, 0, num_atoms * sizeof(cu_real_t));
			// cudaMemset(cu::d_z_hamr_field, 0, num_atoms * sizeof(cu_real_t));

			check_cuda_errors (__FILE__, __LINE__);

			if(::hamr::head_laser_on){

				check_cuda_errors (__FILE__, __LINE__);

				const cu_real_t Tmin = sim::Tmin;
				const cu_real_t Tmax = sim::Tmax;
				const cu_real_t Hx_app = sim::H_vec[0]*sim::H_applied;
				const cu_real_t Hy_app = sim::H_vec[1]*sim::H_applied;
				const cu_real_t Hz_app = sim::H_vec[2]*sim::H_applied;
				// Determine constants
				const cu_real_t H_bounds_x = cu::hamr::d_H_bounds_x;
				const cu_real_t H_bounds_y = cu::hamr::d_H_bounds_y;
				const cu_real_t NPS = cu::hamr::d_NPS;
				const cu_real_t laser_sigma_x2 = cu::hamr::d_laser_sigma_x * cu::hamr::d_laser_sigma_x;
				const cu_real_t laser_sigma_y2 = cu::hamr::d_laser_sigma_y * cu::hamr::d_laser_sigma_y;
				// Update head position - updated in src/hamr/hamr_continuous.cpp
				cu::hamr::d_head_position_x = ::hamr::get_head_position_x();
				cu::hamr::d_head_position_y = ::hamr::get_head_position_y();
				const cu_real_t px = cu::hamr::d_head_position_x;
				const cu_real_t py = cu::hamr::d_head_position_y;

				// Apply thermal field
				apply_local_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::atoms::d_x_coord, cu::atoms::d_y_coord,
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					Tmin, Tmax,
					laser_sigma_x2, laser_sigma_y2,
					px, py,	 
					cu::d_rand_state,
					cu::mp::d_material_params,
					cu::atoms::d_materials,
					num_atoms);

				check_cuda_errors (__FILE__, __LINE__);

				// Apply external field 
				apply_local_external_field_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::atoms::d_x_coord, cu::atoms::d_y_coord,
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					Hx_app, Hy_app, Hz_app,
					px, py,
					H_bounds_x, H_bounds_y, 
					NPS,
					num_atoms);

				check_cuda_errors (__FILE__, __LINE__);
			}
			// Apply global temperature if laser is off
			else{

				check_cuda_errors (__FILE__, __LINE__);

				cu::apply_global_temperature_kernel <<< cu::grid_size, cu::block_size >>> (
					cu::d_x_hamr_field, cu::d_y_hamr_field, cu::d_z_hamr_field,
					global_temperature,
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