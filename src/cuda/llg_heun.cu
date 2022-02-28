//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "exchange_fields.hpp"
#include "internal.hpp"
#include "llg_heun.hpp"

// Namespace aliasing for brevity
#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda{

   //--------------------------------------------------------------------------
   // Function to perform a single heun step
   //--------------------------------------------------------------------------
   void llg_heun(){

      #ifdef CUDA
        // check for cuda initialization, and initialize if necessary
        if (!cu::llg::initialized) cu::llg::__llg_init();
        // perform a single LLG Heun step
        cu::llg::__llg_step();

      #endif

      return;
   }

#ifdef CUDA

	namespace internal {

		//--------------------------------------------------------------
      // device data structures and functions for llg integration
		//--------------------------------------------------------------
      namespace llg{

         // flag to indicate initialization of llg data structures and variables
         bool initialized(false);

         // arrays for storing temporary spin data for integration
         /*cu_real_array_t x_spin_buffer_array(0UL);
         cu_real_array_t y_spin_buffer_array(0UL);
         cu_real_array_t z_spin_buffer_array(0UL);
         cu_real_array_t dS_x_array(0UL);
         cu_real_array_t dS_y_array(0UL);
         cu_real_array_t dS_z_array(0UL);
         */

         // array for storing heun integration prefactors
         //thrust::device_vector<heun_parameters_t> heun_parameters_device(0UL);
         heun_parameters_t *d_heun_params;

         cu_real_t *d_x_spin_buffer;
         cu_real_t *d_y_spin_buffer;
         cu_real_t *d_z_spin_buffer;

         cu_real_t *d_ds_x;
         cu_real_t *d_ds_y;
         cu_real_t *d_ds_z;

         //-----------------------------------------------------------
         // Function to initialize LLG variables on device
         //-----------------------------------------------------------
         void __llg_init (){

            // Reserve space for the buffers
            //cu::llg::x_spin_buffer_array.resize(::atoms::num_atoms);
            //cu::llg::y_spin_buffer_array.resize(::atoms::num_atoms);
            //cu::llg::z_spin_buffer_array.resize(::atoms::num_atoms);
            //cu::llg::dS_x_array.resize(::atoms::num_atoms);
            //cu::llg::dS_y_array.resize(::atoms::num_atoms);
            //cu::llg::dS_z_array.resize(::atoms::num_atoms);

            cudaMalloc((void**)&cu::llg::d_x_spin_buffer, ::atoms::num_atoms * sizeof(cu_real_t));
            cudaMalloc((void**)&cu::llg::d_y_spin_buffer, ::atoms::num_atoms * sizeof(cu_real_t));
            cudaMalloc((void**)&cu::llg::d_z_spin_buffer, ::atoms::num_atoms * sizeof(cu_real_t));

            // Initial copy to the buffer
	    // This is later taken care of inside the corrector step kernel
            cudaMemcpy(cu::llg::d_x_spin_buffer, cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
            cudaMemcpy(cu::llg::d_y_spin_buffer, cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
            cudaMemcpy(cu::llg::d_z_spin_buffer, cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
            
            cudaMalloc((void**)&cu::llg::d_ds_x, ::atoms::num_atoms * sizeof(cu_real_t));
            cudaMalloc((void**)&cu::llg::d_ds_y, ::atoms::num_atoms * sizeof(cu_real_t));
            cudaMalloc((void**)&cu::llg::d_ds_z, ::atoms::num_atoms * sizeof(cu_real_t));

            // Initialize host array for heun parameters
            size_t num_mats = ::mp::num_materials;
            //thrust::host_vector<heun_parameters_t> heun_parameters_host(num_mats);
            std::vector<heun_parameters_t> h_heun_params(num_mats);

            // loop over all materials and determine prefactor constants
            for (size_t i = 0; i < num_mats; i++){

               // temporary variables for readability
               double alpha = ::mp::material[i].alpha;
               double gamma = ::mp::material[i].gamma_rel;

               // save derived variables to host array
               h_heun_params.at(i).prefactor = -gamma / (1.0 + alpha * alpha);
               h_heun_params.at(i).lambda_times_prefactor = -gamma * alpha / (1.0 + alpha * alpha);

               #ifdef CUDA_SPIN_DEBUG
					   // output variables to screen for debugging
                  std::cout << "Heun parameters: "
									 << h_heun_params.at(i).prefactor << " "
									 << h_heun_params.at(i).lambda_times_prefactor << std::endl;
               #endif

            }

            // resize device parameter array to correct size
            //cu::llg::heun_parameters_device.resize(num_mats);
            // copy host array to device
            //thrust::copy(heun_parameters_host.begin(),heun_parameters_host.end(),cu::llg::heun_parameters_device.begin());
            cudaMalloc((void**)&cu::llg::d_heun_params, num_mats * sizeof(heun_parameters_t));
            cudaMemcpy(cu::llg::d_heun_params, h_heun_params.data(), num_mats * sizeof(heun_parameters_t), cudaMemcpyHostToDevice);

            // set flag to indicate data initialization
            initialized=true;

				// check for errors?
				return;

         }

         //---------------------------------------------------------
         // Function to calculate single heun step
         //---------------------------------------------------------
         void __llg_step(){

            // Make a device to device (D2D) copy of initial spin directions
            //thrust::copy (cu::atoms::x_spin_array.begin(),cu::atoms::x_spin_array.end(),cu::llg::x_spin_buffer_array.begin());
            //thrust::copy (cu::atoms::y_spin_array.begin(),cu::atoms::y_spin_array.end(),cu::llg::y_spin_buffer_array.begin());
            //thrust::copy (cu::atoms::z_spin_array.begin(),cu::atoms::z_spin_array.end(),cu::llg::z_spin_buffer_array.begin());

            // D2D copy - is it really necessary?
            //cudaMemcpy(cu::llg::d_x_spin_buffer, cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
            //cudaMemcpy(cu::llg::d_y_spin_buffer, cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);
            //cudaMemcpy(cu::llg::d_z_spin_buffer, cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyDeviceToDevice);

            #ifdef CUDA_SPIN_DEBUG
               // Output first spin position
				   // std::cout << cu::atoms::x_spin_array[0] << " "
					// 			 << cu::atoms::y_spin_array[0] << " "
					// 			 << cu::atoms::z_spin_array[0] << std::endl;
            #endif

            // Cast device pointers for calling kernel function calls
            //cu_real_t * d_x_spin = thrust::raw_pointer_cast(cu::atoms::x_spin_array.data());
            //cu_real_t * d_y_spin = thrust::raw_pointer_cast(cu::atoms::y_spin_array.data());
            //cu_real_t * d_z_spin = thrust::raw_pointer_cast(cu::atoms::z_spin_array.data());

            //int * d_materials = thrust::raw_pointer_cast(cu::atoms::type_array.data());

            //cu::heun_parameters_t * d_heun_params = thrust::raw_pointer_cast (cu::llg::heun_parameters_device.data());

            /*cu_real_t * d_x_spin_field = thrust::raw_pointer_cast(cu::x_total_spin_field_array.data());
            cu_real_t * d_y_spin_field = thrust::raw_pointer_cast(cu::y_total_spin_field_array.data());
            cu_real_t * d_z_spin_field = thrust::raw_pointer_cast(cu::z_total_spin_field_array.data());

            cu_real_t * d_x_external_field = thrust::raw_pointer_cast(cu::x_total_external_field_array.data());
            cu_real_t * d_y_external_field = thrust::raw_pointer_cast(cu::y_total_external_field_array.data());
            cu_real_t * d_z_external_field = thrust::raw_pointer_cast(cu::z_total_external_field_array.data());

            cu_real_t * d_x_spin_buffer = thrust::raw_pointer_cast(cu::llg::x_spin_buffer_array.data());
            cu_real_t * d_y_spin_buffer = thrust::raw_pointer_cast(cu::llg::y_spin_buffer_array.data());
            cu_real_t * d_z_spin_buffer = thrust::raw_pointer_cast(cu::llg::z_spin_buffer_array.data());

            cu_real_t * dS_x_dptr = thrust::raw_pointer_cast(cu::llg::dS_x_array.data());
            cu_real_t * dS_y_dptr = thrust::raw_pointer_cast(cu::llg::dS_y_array.data());
            cu_real_t * dS_z_dptr = thrust::raw_pointer_cast(cu::llg::dS_z_array.data());
            */
            // Calculate spin dependent fields
            cu::update_spin_fields ();

            // Calculate external fields (fixed for integration step)
            cu::update_external_fields ();

            // Check for cuda errors in file, line
            check_cuda_errors (__FILE__, __LINE__);

            // Invoke kernel for heun predictor step
            cu::llg::llg_heun_predictor_step <<< cu::grid_size, cu::block_size >>> (
               cu::atoms::d_materials, cu::llg::d_heun_params,
					cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
               cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
               cu::llg::d_ds_x, cu::llg::d_ds_y, cu::llg::d_ds_z,
               ::mp::dt, ::atoms::num_atoms);

            #ifdef CUDA_SPIN_DEBUG
               // std::cout << cu::atoms::x_spin_array[0] << " "
               //           << cu::atoms::y_spin_array[0] << " "
               //           << cu::atoms::z_spin_array[0] << " "
               //           << cu::z_total_spin_field_array[0] << " "
               //           << cu::z_total_external_field_array[0] << " "
               //           << std::endl;
            #endif

            check_cuda_errors (__FILE__, __LINE__);

            // Recalculate spin dependent fields
            cu::update_spin_fields ();

            check_cuda_errors (__FILE__, __LINE__);

            // Invoke kernel for heun corrector step
            cu::llg::llg_heun_corrector_step <<< cu::grid_size, cu::block_size >>> (
               cu::atoms::d_materials, cu::llg::d_heun_params,
					cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
               cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
               cu::d_x_external_field, cu::d_y_external_field, cu::d_z_external_field,
               cu::llg::d_x_spin_buffer, cu::llg::d_y_spin_buffer, cu::llg::d_z_spin_buffer,
               cu::llg::d_ds_x, cu::llg::d_ds_y, cu::llg::d_ds_z,
               ::mp::dt, ::atoms::num_atoms);

            check_cuda_errors (__FILE__, __LINE__);

				return;

         }

         //-----------------------------------------------------------------------------
         // CUDA Kernel to calculate the predictor step of the Landau-Lifshitz-Gilbert
         // (LLG) equation using the Heun scheme
         //-----------------------------------------------------------------------------
         __global__ void llg_heun_predictor_step (
               int * material_id,
               cu::heun_parameters_t * heun_parameters,
               cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin, // contains initial spin
               cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
               cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
               cu_real_t * dSx, cu_real_t * dSy, cu_real_t * dSz,
               cu_real_t dt, size_t num_atoms){

            // Loop over blocks for large systems > ~100k spins
            for ( size_t atom = blockIdx.x * blockDim.x + threadIdx.x;
                  atom < num_atoms;
                  atom += blockDim.x * gridDim.x)
            {

               // Determine material id for atom
               size_t mid = material_id[atom];

               // prestore prefactors into registers
               const cu_real_t prefactor = heun_parameters[mid].prefactor;
               const cu_real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

               // load spin direction to registers for later multiple reuse
               cu_real_t sx = x_spin[atom];
               cu_real_t sy = y_spin[atom];
               cu_real_t sz = z_spin[atom];

               // calculate total field
               cu_real_t H_x = x_sp_field[atom] + x_ext_field[atom];
               cu_real_t H_y = y_sp_field[atom] + y_ext_field[atom];
               cu_real_t H_z = z_sp_field[atom] + z_ext_field[atom];

               // calculate s x h
               cu_real_t sxh_x = sy * H_z - sz * H_y;
               cu_real_t sxh_y = sz * H_x - sx * H_z;
               cu_real_t sxh_z = sx * H_y - sy * H_x;

               // calculate delta S
               cu_real_t Ds_x = prefactor * sxh_x + lambdatpr * (sy * sxh_z - sz * sxh_y); // sxsxh_x;
               cu_real_t Ds_y = prefactor * sxh_y + lambdatpr * (sz * sxh_x - sx * sxh_z); // sxsxh_y;
               cu_real_t Ds_z = prefactor * sxh_z + lambdatpr * (sx * sxh_y - sy * sxh_x); // sxsxh_z;

               // store initial change in S for second step
               dSx[atom] = Ds_x;
               dSy[atom] = Ds_y;
               dSz[atom] = Ds_z;

               // Calculate intermediate spin position
               cu_real_t new_spin_x = sx + Ds_x * dt;
               cu_real_t new_spin_y = sy + Ds_y * dt;
               cu_real_t new_spin_z = sz + Ds_z * dt;

               // calculate spin length for renormalization
               #ifdef CUDA_DP
                  double mod_s = 1.0 / __dsqrt_rn(
               #else
                  // cuda intrinsic precise reciprocal square root
                  float mod_s = __frsqrt_rn(
               #endif
                  new_spin_x * new_spin_x +
                  new_spin_y * new_spin_y +
                  new_spin_z * new_spin_z
               );

               // Store normalized intermediate spin direction
               x_spin[atom] = new_spin_x * mod_s;
               y_spin[atom] = new_spin_y * mod_s;
               z_spin[atom] = new_spin_z * mod_s;

            }

         }

         //-----------------------------------------------------------------------------
			// CUDA Kernel to calculate the corrector step of the Landau-Lifshitz-Gilbert
			// (LLG) equation using the Heun scheme
			//-----------------------------------------------------------------------------
         __global__ void llg_heun_corrector_step (
               int * material_id,
               cu::heun_parameters_t * heun_parameters,
               cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin, // spin from predictor step
               cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
               cu_real_t * x_ext_field, cu_real_t * y_ext_field, cu_real_t * z_ext_field,
               cu_real_t * x_spin_buffer, cu_real_t * y_spin_buffer, cu_real_t * z_spin_buffer, // initial spin
               cu_real_t * dSx, cu_real_t * dSy, cu_real_t * dSz,
               cu_real_t dt, size_t num_atoms
               )
         {

            for ( size_t atom = blockIdx.x * blockDim.x + threadIdx.x;
                  atom < num_atoms;
                  atom += blockDim.x * gridDim.x)
            {

               size_t mid = material_id[atom];

               cu_real_t prefactor = heun_parameters[mid].prefactor;
               cu_real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

               // copy the predictor spins to registers
               cu_real_t spin_x = x_spin[atom];
               cu_real_t spin_y = y_spin[atom];
               cu_real_t spin_z = z_spin[atom];

               //the field
               cu_real_t H_x = x_sp_field[atom] + x_ext_field[atom];
               cu_real_t H_y = y_sp_field[atom] + y_ext_field[atom];
               cu_real_t H_z = z_sp_field[atom] + z_ext_field[atom];

               // calculate components
               cu_real_t SxH_x = spin_y * H_z - spin_z * H_y;
               cu_real_t SxH_y = spin_z * H_x - spin_x * H_z;
               cu_real_t SxH_z = spin_x * H_y - spin_y * H_x;

               cu_real_t SxSxH_x = spin_y * SxH_z - spin_z * SxH_y;
               cu_real_t SxSxH_y = spin_z * SxH_x - spin_x * SxH_z;
               cu_real_t SxSxH_z = spin_x * SxH_y - spin_y * SxH_x;

               // Calculate dS
               cu_real_t DS_prime_x = prefactor * SxH_x + lambdatpr * SxSxH_x;
               cu_real_t DS_prime_y = prefactor * SxH_y + lambdatpr * SxSxH_y;
               cu_real_t DS_prime_z = prefactor * SxH_z + lambdatpr * SxSxH_z;

               // Calculate heun step
               cu_real_t S_x = x_spin_buffer[atom] + 0.5 * (dSx[atom] + DS_prime_x) * dt;
               cu_real_t S_y = y_spin_buffer[atom] + 0.5 * (dSy[atom] + DS_prime_y) * dt;
               cu_real_t S_z = z_spin_buffer[atom] + 0.5 * (dSz[atom] + DS_prime_z) * dt;

               #ifdef CUDA_DP
                  double mods = 1.0 / __dsqrt_rn(S_x*S_x + S_y*S_y + S_z*S_z);
               #else
                  // cuda intrinsic precise reciprocal square root
                  float mods = __frsqrt_rn(S_x*S_x + S_y*S_y + S_z*S_z);
               #endif

	       cu_real_t sp_x_new = mods * S_x;
	       cu_real_t sp_y_new = mods * S_y;
	       cu_real_t sp_z_new = mods * S_z;

               // save final spin direction to spin array
               x_spin[atom] = sp_x_new;
               y_spin[atom] = sp_y_new;
               z_spin[atom] = sp_z_new;

	       // Update the spin buffer
               x_spin_buffer[atom] = sp_x_new;
               y_spin_buffer[atom] = sp_y_new;
               z_spin_buffer[atom] = sp_z_new;

            }
         }
      } /* llg */
   }

#endif

} // end of namespace cuda
