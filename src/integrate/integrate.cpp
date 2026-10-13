//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
//   Serial and parallel integrator dispatch, moved from src/simulate/sim.cpp
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cstdint>
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "environment.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "montecarlo.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "../cells/internal.hpp"
#include "../micromagnetic/internal.hpp"

// integrate module headers
#include "integrate.hpp"
#include "internal.hpp"

namespace integrate{

//------------------------------------------------------------------------------
// Wrapper function to call serial or parallel integrators
//------------------------------------------------------------------------------
int integrate(uint64_t n_steps){

	// Check for calling of function
	if(err::check==true) std::cout << "integrate::integrate has been called" << std::endl;

	// Call serial or parallell depending at compile time
	#ifdef MPICF
		internal::integrate_mpi(n_steps);
	#else
		internal::integrate_serial(n_steps);
	#endif

	// return
	return EXIT_SUCCESS;
}

namespace internal{

/// @brief Wrapper function to call serial integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls serial integrators based on integrate::integrator
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.2
/// @date    23/04/2021
///
/// @internal
///	Created:		05/02/2011
///	Revision:	    23/04/2021
///=====================================================================================
///
void integrate_serial(uint64_t n_steps){

   // Check for calling of function
   if(err::check==true) std::cout << "integrate::internal::integrate_serial has been called" << std::endl;

	// if simulation is micromagnetic
   if (micromagnetic::discretisation_type >0 ) multiscale_simulation_steps(n_steps);


   //else simulation is atomistic
   else{

   // Case statement to call integrator
   switch(integrate::integrator){

      case integrator_t::llg_heun: // LLG Heun
         for(uint64_t ti=0;ti<n_steps;ti++){
            // Optionally select GPU accelerated version
            if(gpu::acceleration) gpu::llg_heun();
            // Otherwise use CPU version
            else internal::LLG_Heun();
				if (environment::enabled && (sim::time)%environment::num_atomic_steps_env ==0){
					environment::LLB(sim::temperature, sim::H_applied,sim::H_vec[0],sim::H_vec[1],sim::H_vec[2],mp::dt);
				}
            // Increment time
            internal::increment_time();
         }
         break;

		case integrator_t::monte_carlo: // Montecarlo
			for(uint64_t ti=0;ti<n_steps;ti++){

                // Optionally select GPU accelerated version
                if(gpu::acceleration) gpu::mc_step();

                else montecarlo::mc_step(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::num_atoms, atoms::type_array);

				// increment time
				internal::increment_time();
			}
			break;

      case integrator_t::llg_midpoint: // LLG Midpoint
         for(uint64_t ti=0;ti<n_steps;ti++){
            internal::LLG_Midpoint();
            // increment time
            internal::increment_time();
         }
         break;

		case integrator_t::cmc: // Constrained Monte Carlo
			for(uint64_t ti=0;ti<n_steps;ti++){
				montecarlo::cmc_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::hybrid_cmc: // Hybrid Constrained Monte Carlo
			for(uint64_t ti=0;ti<n_steps;ti++){
				montecarlo::cmc_mc_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_quantum: // LLG quantum noise
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::llg_quantum_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_rkmk2: // LLG Runge-Kutta-Munthe-Kaas 2nd order
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::llg_rkmk2_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_rkmk4: // LLG Runge-Kutta-Munthe-Kaas 4th order
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::llg_rkmk4_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_dm: // LLG Depondt-Mertens
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::llg_dm_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf: // LSF
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::lsf_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf_mc: // LSF Monte Carlo
			for(uint64_t ti=0;ti<n_steps;ti++){
				montecarlo::lsf_mc_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf_rk4: // LSF-RK4
			for(uint64_t ti=0;ti<n_steps;ti++){
				internal::lsf_rk4_step();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::suzuki_trotter: // spin-lattice Dynamics
			for(uint64_t ti=0;ti<n_steps;ti++){
				sld::suzuki_trotter();
			   internal::increment_time();
			}
			break;

		default:{
			std::cerr << "Unknown integrator type "<< static_cast<int>(integrate::integrator) << " requested, exiting" << std::endl;
         err::vexit();
		}
	}
}

   return;
}

/// @brief Wrapper function to call MPI parallel integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls parallel integrators based on integrate::integrator
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		07/03/2011
///	Revision:	  ---
///=====================================================================================
///
int integrate_mpi(uint64_t n_steps){

	// Check for calling of function
	if(err::check==true) std::cout << "integrate::internal::integrate_mpi has been called" << std::endl;


	if (micromagnetic::discretisation_type >0 ) multiscale_simulation_steps(n_steps);

		//else simulation is atomistic
	else{

	// Case statement to call integrator
	switch(integrate::integrator){
		case integrator_t::llg_heun: // LLG Heun
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
				// Select CUDA version if supported
				#ifdef CUDA
					//sim::LLG_Heun_cuda_mpi();
				#else
					internal::LLG_Heun_mpi();
					//calcualte the field from the environment
					if (environment::enabled &&  (sim::time)%environment::num_atomic_steps_env ==0)
						environment::LLB(sim::temperature,
												sim::H_applied,
												sim::H_vec[0],
												sim::H_vec[1],
												sim::H_vec[2],
												mp::dt);
				#endif
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::monte_carlo: // Montecarlo

			for(uint64_t ti=0;ti<n_steps;ti++){
				#ifdef MPICF
               if(montecarlo::mc_parallel_initialized == false) {
                  montecarlo::mc_parallel_init(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
                                               vmpi::min_dimensions, vmpi::max_dimensions);
               }
               montecarlo::mc_step_parallel(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                            atoms::type_array);
            #endif

				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_midpoint: // LLG Midpoint
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
			// Select CUDA version if supported
				#ifdef CUDA
					//sim::LLG_Midpoint_cuda_mpi();
				#else
					internal::LLG_Midpoint_mpi();
				#endif
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_rkmk2: // LLG Runge-Kutta-Munthe-Kaas 2nd order
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
				internal::llg_rkmk2_mpi_step();
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_rkmk4: // LLG Runge-Kutta-Munthe-Kaas 4th order
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
				internal::llg_rkmk4_mpi_step();
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::llg_dm: // LLG Depondt-Mertens
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
				internal::llg_dm_mpi_step();
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::cmc: // Constrained Monte Carlo
			for(uint64_t ti=0;ti<n_steps;ti++){
				terminaltextcolor(RED);
				std::cerr << "Error - Constrained Monte Carlo Integrator unavailable for parallel execution" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf: // LSF
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
			// Select CUDA version if supported
				#ifdef CUDA
					//sim::LSF_cuda();
				#else
					internal::LSF_mpi();
				#endif
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf_mc: // LSF-Montecarlo
			for(uint64_t ti=0;ti<n_steps;ti++){
				#ifdef MPICF
					if(montecarlo::lsf_mc_parallel_initialized == false) {
						montecarlo::lsf_mc_parallel_init(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
																	vmpi::min_dimensions, vmpi::max_dimensions);
					}
					montecarlo::lsf_mc_step_parallel(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
																atoms::type_array);
					#endif

				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::lsf_rk4: // LSF-RK4
			for(uint64_t ti=0;ti<n_steps;ti++){
			#ifdef MPICF
			// Select CUDA version if supported
				#ifdef CUDA
					//sim::LSF_RK4_cuda();
				#else
					internal::LSF_RK4_mpi();
				#endif
			#endif
				// increment time
				internal::increment_time();
			}
			break;

		case integrator_t::suzuki_trotter: // Suzuki Trotter decomposition
 			for(uint64_t ti=0;ti<n_steps;ti++){
 				#ifdef MPICF

                sld::suzuki_trotter_step_parallel(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                             atoms::type_array);
             #endif
 				// increment time
 				internal::increment_time();
 			}
 			break;

		default:{
			terminaltextcolor(RED);
			std::cerr << "Unknown integrator type "<< static_cast<int>(integrate::integrator) << " requested, exiting" << std::endl;
			terminaltextcolor(WHITE);
			exit (EXIT_FAILURE);
			}
		}
	}

	return EXIT_SUCCESS;
}

void multiscale_simulation_steps(int n_steps){

	if (environment::enabled && (sim::time)%environment::num_atomic_steps_env ==0)
		environment::LLB(sim::temperature,
						   	sim::H_applied,
						   	sim::H_vec[0],
						   	sim::H_vec[1],
						   	sim::H_vec[2],
						   	mp::dt);

   for(int ti=0;ti<n_steps;ti++){
      //calcaulte the field from the environment
      // if (environment::enabled && (sim::time)%environment::num_atomic_steps_env ==0)
		//	environment::LLB(sim::temperature,
      //    sim::H_applied,
      //    sim::H_vec[0],
      //    sim::H_vec[1],
      //    sim::H_vec[2],
      //    mp::dt);

   //if  there are micromagnetic cells run a micromagnetic step
   if (micromagnetic::number_of_micromagnetic_cells > 0 &&  (sim::time)% micromagnetic::num_atomic_steps_mm == 0) {


      //if LLb run an LLG steps
      if (micromagnetic::integrator == 0) micromagnetic::LLG(cells::local_cell_array,
         n_steps,
         cells::num_cells,
         cells::num_local_cells,
         sim::temperature,
         cells::mag_array_x,
         cells::mag_array_y,
         cells::mag_array_z,
         sim::H_vec[0],
         sim::H_vec[1],
         sim::H_vec[2],
         sim::H_applied,
         mp::dt,
         cells::volume_array);

         //if LLB run an LLB step

         else micromagnetic::LLB(cells::local_cell_array,
            n_steps,
            cells::num_cells,
            cells::num_local_cells,
            sim::temperature,
            cells::mag_array_x,
            cells::mag_array_y,
            cells::mag_array_z,
            sim::H_vec[0],
            sim::H_vec[1],
            sim::H_vec[2],
            sim::H_applied,
            mp::dt,
            cells::volume_array);
         }
         //run an atomistic step if there are atomistic atoms
         if (micromagnetic::number_of_atomistic_atoms > 0) micromagnetic::atomistic_LLG_Heun();

         //incremenet time
         internal::increment_time();
      }

   }

} // end of internal namespace
} // end of integrate namespace
