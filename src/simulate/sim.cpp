///
/// @file
/// @brief Contains the sim namespace and wrapper functions for system integration
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    09/03/2011
/// @internal
///	Created:		09/03/2011
///	Revision:	  ---
///=====================================================================================
///

// Vampire Header files
#include "atoms.hpp"
#include "program.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <iostream>

namespace sim{
	std::ofstream mag_file;
	int time=0;
	int total_time=1;
	int loop_time=1;
	int partial_time=1;
	int equilibration_time=0;
	int runs=1; // for certain repetitions in programs
	
	bool ext_demag=false;
	
	double Tmax=300.0;
	double Tmin=0.0;
	double Teq=300.0;
	double temperature=300.0;
	double delta_temperature=10.0;
	double H_applied=0.0;
	double H_vec[3]={0.0,0.0,1.0};
	double Hmin=-1.0; // T
	double Hmax=+1.0; // T
	double Hinc= 0.1; // T
	double Heq=0.0;
	double demag_factor[3]={0.0,0.0,0.0};
	double constraint_phi=0.0; // Constrained minimisation vector (azimuthal) [degrees]
	double constraint_theta=0.0; // Constrained minimisation vector (rotational) [degrees]
	double head_position[2]={0.0,cs::system_dimensions[1]*0.5}; // A
	double head_speed=30.0; // nm/ns
	bool   head_laser_on=false;

	double cooling_time=100.0e-12; //seconds
	int cooling_function_flag=0; // 0 = exp, 1 = gaussian
	double pump_power=2.4e22;
	double pump_time=20.0e-15; 
	
	int system_simulation_flags;
	int hamiltonian_simulation_flags[10];
	int integrator=0; // 0 = LLG Heun; 1= MC; 2 = LLG Midpoint; 3 = CMC 
	int program=0; 
	
	bool surface_anisotropy=false; // flag to enable surface anisotropy
	bool identify_surface_atoms=true; // flag to idenify surface atoms in config coordinate file
	unsigned int surface_anisotropy_threshold=6; // threshold for surface atoms

	
	// Local function declarations
	int integrate_serial(int);
	int integrate_mpi(int);
	
/// @brief Function to increment time counter and associted variables
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		02/10/2008
///	Revision:	1.1 09/03/2011
///=====================================================================================
///
	void increment_time(){
		
		sim::time++;
		sim::head_position[0]+=sim::head_speed*mp::dt_SI*1.0e10;
		
	}
	
/// @brief Function to run one a single program
///
/// @callgraph
/// @callergraph
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		02/10/2008
///	Revision:	1.1 09/03/2011
///=====================================================================================
///
int run(){
	// Check for calling of function
	if(err::check==true) std::cout << "sim::run has been called" << std::endl;

	// For MPI version, calculate initialisation time
	if(vmpi::my_rank==0){
		#ifdef MPICF
			std::cout << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
			vmpi::start_time=MPI_Wtime(); // reset timer
		#endif
		std::cout << "Starting Simulation with Program ";
	}
	
	//program::timestep_scaling();
	//return(EXIT_SUCCESS);
	
	// Initialise random number generator
	mtrandom::grnd.seed(mtrandom::integration_seed+vmpi::my_rank);

	
	// Select program to run
	switch(sim::program){
		case 0:
			if(vmpi::my_rank==0) std::cout << "Benchmark..." << std::endl; 
			program::bmark();
			break;
		
		case 1:
			if(vmpi::my_rank==0) std::cout << "Time-Series..." << std::endl; 
			program::time_series();
			break;
		
		case 2: 
			if(vmpi::my_rank==0) std::cout << "Hysteresis-Loop..." << std::endl; 
			program::hysteresis();
			break;
			
		case 3: 
			if(vmpi::my_rank==0) std::cout << "Static-Hysteresis-Loop..." << std::endl; 
			program::static_hysteresis();
			break;
			
		case 4:
			if(vmpi::my_rank==0) std::cout << "Curie-Temperature..." << std::endl; 
			program::curie_temperature();
			break;
			
		case 5:
			if(vmpi::my_rank==0) std::cout << "Field-Cool..." << std::endl; 
			program::field_cool();
			break;

		case 6:
			if(vmpi::my_rank==0) std::cout << "Two-Temperature-Pulse..." << std::endl; 
			program::two_temperature_pulse();
			break;
			
		case 7:
			if(vmpi::my_rank==0) std::cout << "HAMR-Simulation..." << std::endl; 
			program::hamr();
			break;
		
		case 50:
			if(vmpi::my_rank==0) std::cout << "Diagnostic-Boltzmann..." << std::endl; 
			program::boltzmann_dist();
			break;
		
		default:{
			std::cerr << "Unknown Internal Program ID "<< sim::program << " requested, exiting" << std::endl;
			exit (EXIT_FAILURE);
			}
	}

	// output Monte Carlo Statistics if applicable
	if(sim::integrator==3){
		std::cout << "Constrained Monte Carlo Statistics:" << std::endl;
		std::cout << "\tTotal moves: " << cmc::mc_total << std::endl;
		std::cout << "\t" << (cmc::mc_success/cmc::mc_total)*100.0 << "% Accepted" << std::endl;
		std::cout << "\t" << (cmc::energy_reject/cmc::mc_total)*100.0 << "% Rejected (Energy)" << std::endl;
		std::cout << "\t" << (cmc::sphere_reject/cmc::mc_total)*100.0 << "% Rejected (Sphere)" << std::endl;
		
	}

	//program::LLB_Boltzmann();

	return EXIT_SUCCESS;
}

/// @brief Wrapper function to call integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls serial or parallel integrator 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int integrate(int n_steps){
	
	// Check for calling of function
	if(err::check==true) std::cout << "sim::integrate has been called" << std::endl;
	
	// Call serial or parallell depending at compile time
	#ifdef MPICF
		sim::integrate_mpi(n_steps);
	#else 
		sim::integrate_serial(n_steps);
	#endif
	
	// return
	return EXIT_SUCCESS;
}

/// @brief Wrapper function to call serial integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls serial integrators based on sim::integrator 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int integrate_serial(int n_steps){
	
	// Check for calling of function
	if(err::check==true) std::cout << "sim::integrate_serial has been called" << std::endl;
	
	// Case statement to call integrator
	switch(sim::integrator){
		case 0: // LLG Heun
			for(int ti=0;ti<n_steps;ti++){
				// Select CUDA version if supported
				#ifdef CUDA
					sim::LLG_Heun_cuda();
				#else
					sim::LLG_Heun();
				#endif
				// increment time
				increment_time();
			}
			break;
		
		case 1: // Montecarlo
			for(int ti=0;ti<n_steps;ti++){
				sim::MonteCarlo();
				// increment time
				increment_time();
			}
			break;
		
		case 2: // LLG Midpoint
			for(int ti=0;ti<n_steps;ti++){
				// Select CUDA version if supported
				#ifdef CUDA
					sim::LLG_Midpoint_cuda();
				#else
					sim::LLG_Midpoint();
				#endif
				// increment time
				increment_time();
			}
			break;
			
		case 3: // Constrained Monte Carlo
			for(int ti=0;ti<n_steps;ti++){
				sim::ConstrainedMonteCarlo();
				// increment time
				increment_time();
			}
			break;

		default:{
			std::cerr << "Unknown integrator type "<< sim::integrator << " requested, exiting" << std::endl;
			exit (EXIT_FAILURE);
			}
	}
	
	return EXIT_SUCCESS;
}

/// @brief Wrapper function to call MPI parallel integrators
///
/// @callgraph
/// @callergraph
///
/// @details Calls parallel integrators based on sim::integrator 
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
int integrate_mpi(int n_steps){
	
	// Check for calling of function
	if(err::check==true) std::cout << "sim::integrate_mpi has been called" << std::endl;
	
	// Case statement to call integrator
	switch(sim::integrator){
		case 0: // LLG Heun
			for(int ti=0;ti<n_steps;ti++){
			#ifdef MPICF
				// Select CUDA version if supported
				#ifdef CUDA
					//sim::LLG_Heun_cuda_mpi();
				#else
					sim::LLG_Heun_mpi();
				#endif
			#endif
				// increment time
				increment_time();
			}
			break;
		
		case 1: // Montecarlo
			for(int ti=0;ti<n_steps;ti++){
				std::cerr << "Error - Monte Carlo Integrator unavailable for parallel execution" << std::endl;
				err::vexit();
				// increment time
				increment_time();
			}
			break;
		
		case 2: // LLG Midpoint
			for(int ti=0;ti<n_steps;ti++){
			#ifdef MPICF
			// Select CUDA version if supported
				#ifdef CUDA
					//sim::LLG_Midpoint_cuda_mpi();
				#else
					sim::LLG_Midpoint_mpi();
				#endif
			#endif
				// increment time
				increment_time();
			}
			break;
			
		case 3: // Constrained Monte Carlo
			for(int ti=0;ti<n_steps;ti++){
				std::cerr << "Error - Constrained Monte Carlo Integrator unavailable for parallel execution" << std::endl;
				err::vexit();
				// increment time
				increment_time();
			}
			break;
			
		default:{
			std::cerr << "Unknown integrator type "<< sim::integrator << " requested, exiting" << std::endl;
			exit (EXIT_FAILURE);
			}
	}
	
	return EXIT_SUCCESS;
}

} // Namespace sim


