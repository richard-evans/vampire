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

	
	int run(){
		
	if(vmpi::my_rank==0){
		#ifdef MPICF
			std::cout << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
			vmpi::start_time=MPI_Wtime(); // reset timer
		#endif
		std::cout << "Starting Simulation..." << std::endl;
	}
    //program::curie_temperature(true);
	//program::hamr_run();
	//program::static_hysteresis();
	//program::two_temperature_pulse();
	program::bmark();
	//program::LLB_Boltzmann();
	//program::hysteresis();
	return EXIT_SUCCESS;
}
	// derived variables
	
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


