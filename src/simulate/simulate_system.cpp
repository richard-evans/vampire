//====================================================================================================
//
//       				                    simulate_system
//
//  			 Subroutine to simulate an atomistic system with predefined integration scheme
//				 simulation time, temperature etc
//	 
//								Version 1.0 R Evans 02/10/2008
//
//==================================================================================================== 
#include "atoms.hpp"
#include "program.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include <iostream>

	int set_LLG();
	int set_demag();

namespace sim{
	std::ofstream mag_file;
	int time;
	int total_time=1;
	int loop_time=1;
	int partial_time=1;
	int equilibration_time=0;
	
	double Tmax=300;
	double temperature;
        double delta_temperature;
	double H_applied=0.0;
	double H_vec[3]={0.0,0.0,1.0};
	double Hmin=-1.0; // T
	double Hmax=+1.0; // T
	double Hinc= 0.1; // T	

	int system_simulation_flags;
	int hamiltonian_simulation_flags[10];
	int integrator=0;
	
	// Local function declarations
	int integrate_serial(int);
	int integrate_mpi(int);
	
	int run(){
		
	if(vmpi::my_rank==0){
		#ifdef MPICF
			std::cout << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
			vmpi::start_time=MPI_Wtime(); // reset timer
		#endif
		std::cout << "Starting Simulation..." << std::endl;
	}
    program::curie_temperature(true);
	//program::hamr_run();
	//program::static_hysteresis();
	//program::two_temperature_pulse();
	//program::bmark();
	//program::LLB_Boltzmann();
	//program::hysteresis();
	return EXIT_SUCCESS;
}
	// derived variables

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
		//sim::integrate_mpi(n_steps, istart, iend);
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
			}
			break;
		
		case 1: // Montecarlo
			for(int ti=0;ti<n_steps;ti++){
				sim::MonteCarlo();
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
			}
			break;
			
		case 3: // Constrained Monte Carlo
			for(int ti=0;ti<n_steps;ti++){
				sim::ConstrainedMonteCarlo();
			}
			break;
			
		case 4: // Grain Growth Method
			//grain_growth(cs_num_atoms,cs_coord_array,particle_include_array,cs_atom_type_array);
			std::cerr << "Grain growth not yet implemented, exiting" << std::endl;
			exit(0);
			break;
			
		default:{
			std::cerr << "Unknown integrator type "<< sim::integrator << " requested, exiting" << std::endl;
			exit (EXIT_FAILURE);
			}
	}
	
	// return
	return EXIT_SUCCESS;
}




int initialise(){

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "initialise_system has been called" << std::endl;}

    for(int atom=0;atom<=atoms::num_atoms-1;atom++){
		atoms::x_spin_array[atom] = 0.0;
		atoms::y_spin_array[atom] = 0.1;
		atoms::z_spin_array[atom] = 0.9;
	}
	//std::cout.setf(std::ios::fixed,std::ios::floatfield);
	
  	//sim::mag_file.open ("M_vs_T.txt");
      
	

	set_LLG();
	
	return 0;

}

} // Namespace sim


