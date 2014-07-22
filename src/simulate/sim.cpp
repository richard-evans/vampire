//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
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
#include "demag.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <iostream>

namespace sim{
	std::ofstream mag_file;
	uint64_t time=0;
	uint64_t total_time=10000;
	uint64_t loop_time=10000;
	int partial_time=1000;
	uint64_t equilibration_time=0;
	int runs=1; /// for certain repetitions in programs
	
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
	double applied_field_angle_phi=0.0;
	double applied_field_angle_theta=0.0;
	bool applied_field_set_by_angle=false;
	
	double demag_factor[3]={0.0,0.0,0.0};
	double head_position[2]={0.0,cs::system_dimensions[1]*0.5}; // A
	double head_speed=30.0; /// nm/ns
	bool   head_laser_on=false;
	bool   constraint_rotation=false; /// enables rotation of spins to new constraint direction
	bool   constraint_phi_changed=false; /// flag to note change in phi
	double constraint_phi=0.0; /// Constrained minimisation vector (azimuthal) [degrees]
	double constraint_phi_min=0.0; /// loop angle min [degrees]
	double constraint_phi_max=0.0; /// loop angle max [degrees]
	double constraint_phi_delta=5.0; /// loop angle delta [degrees]

	bool   constraint_theta_changed=false;
	double constraint_theta=0.0; /// Constrained minimisation vector (rotational) [degrees]
	double constraint_theta_min=0.0; /// loop angle min [degrees]
	double constraint_theta_max=0.0; // loop angle max [degrees]
	double constraint_theta_delta=5.0; /// loop angle delta [degrees]
	
	// LaGrange multiplier variables
	double lagrange_lambda_x=0.0;
   double lagrange_lambda_y=0.0;
   double lagrange_lambda_z=0.0;
   double lagrange_m=1.0;
   double lagrange_N=1000.0;
   bool   lagrange_multiplier=false;

	double cooling_time=100.0e-12; ///seconds
	int cooling_function_flag=0; /// 0 = exp, 1 = gaussian
	pump_functions_t pump_function=two_temperature;
	double pump_power=2.4e22;
	double pump_time=20.0e-15; 
	double double_pump_power=2.2e22;
	double double_pump_Tmax=500.0;
	double double_pump_time=10.0e-15; 
	double double_pump_delay=10.0e-12;
	double HeatSinkCouplingConstant=0.0; ///1.1e12 ~ sensible value
	double TTCe = 5.6E04; ///electron specific heat
	double TTCl = 3.0E06; ///phonon specific heat
	double TTG = 17.0E17 ;///electron coupling constant
	double TTTe = 0.0; /// electron temperature
	double TTTp = 0.0; /// phonon temperature
  
   double mc_delta_angle=0.1; /// Tuned angle for Monte Carlo trial move
   mc_algorithms mc_algorithm=hinzke_nowak;
  
	int system_simulation_flags;
	int hamiltonian_simulation_flags[10];
	int integrator=0; /// 0 = LLG Heun; 1= MC; 2 = LLG Midpoint; 3 = CMC 
	int program=0; 
	int AnisotropyType=2; /// Controls scalar (0) or tensor(1) anisotropy (off(2))
	
	bool surface_anisotropy=false; /// flag to enable surface anisotropy
	bool identify_surface_atoms=false; /// flag to idenify surface atoms in config coordinate file
	unsigned int surface_anisotropy_threshold=123456789; /// global threshold for surface atoms
	bool NativeSurfaceAnisotropyThreshold=false; /// enables site-dependent surface threshold
	
	// Anisotropy control booleans
	bool UniaxialScalarAnisotropy=false; /// Enables scalar uniaxial anisotropy
	bool TensorAnisotropy=false; /// Overrides scalar uniaxial anisotropy
	bool second_order_uniaxial_anisotropy=false; /// Enables second order uniaxial anisotropy
   bool sixth_order_uniaxial_anisotropy=false; /// Enables sixth order uniaxial anisotropy
	bool CubicScalarAnisotropy=false; /// Enables scalar cubic anisotropy
   bool lattice_anisotropy_flag=false; /// Enables lattice anisotropy

	bool local_temperature=false; /// flag to enable material specific temperature
	bool local_applied_field=false; /// flag to enable material specific applied field
	bool local_fmr_field=false; /// flag to enable material specific fmr field

   // Checkpoint flags and variables
   bool load_checkpoint_flag=false; // Load spin configurations
   bool load_checkpoint_continue_flag=true; // Continue simulation from checkpoint time
   bool save_checkpoint_flag=false; // Save checkpoint
   bool save_checkpoint_continuous_flag=false; // save checkpoints during simulations
   int save_checkpoint_rate=1; // Default increment between checkpoints

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
		if(sim::hamiltonian_simulation_flags[4]==1) demag::update();
		if(sim::lagrange_multiplier) update_lagrange_lambda();
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
			zlog << zTs() << "Time for initialisation: " << MPI_Wtime()-vmpi::start_time << std::endl;
			vmpi::start_time=MPI_Wtime(); // reset timer
		#endif
		std::cout << "Starting Simulation with Program ";
		zlog << zTs() << "Starting Simulation with Program ";
	}

	// Now set initial compute time
	#ifdef MPICF
	vmpi::ComputeTime=MPI_Wtime();
	vmpi::WaitTime=MPI_Wtime();
	vmpi::TotalComputeTime=0.0;
	vmpi::TotalWaitTime=0.0;
	#endif

	// Initialise random number generator
	mtrandom::grnd.seed(mtrandom::integration_seed+vmpi::my_rank);

   // Seeds with single bit differences are not ideal and may be correlated for first few values - warming up integrator
   for(int i=0; i<1000; ++i) mtrandom::grnd();

   // Check for load spin configurations from checkpoint
   if(sim::load_checkpoint_flag) load_checkpoint();

	// Select program to run
	switch(sim::program){
		case 0:
			if(vmpi::my_rank==0){
				std::cout << "Benchmark..." << std::endl;
				zlog << "Benchmark..." << std::endl;
			}
			program::bmark();
			break;
		
		case 1:
			if(vmpi::my_rank==0){
				std::cout << "Time-Series..." << std::endl;
				zlog << "Time-Series..." << std::endl;
			}
			program::time_series();
			break;
		
		case 2: 
			if(vmpi::my_rank==0){
				std::cout << "Hysteresis-Loop..." << std::endl;
				zlog << "Hysteresis-Loop..." << std::endl;
			}
			program::hysteresis();
			break;
			
		case 3: 
			if(vmpi::my_rank==0){
				std::cout << "Static-Hysteresis-Loop..." << std::endl;
				zlog << "Static-Hysteresis-Loop..." << std::endl;
			}
			program::static_hysteresis();
			break;
			
		case 4:
			if(vmpi::my_rank==0){
				std::cout << "Curie-Temperature..." << std::endl;
				zlog << "Curie-Temperature..." << std::endl;
			}
			program::curie_temperature();
			break;
			
		case 5:
			if(vmpi::my_rank==0){
				std::cout << "Field-Cool..." << std::endl;
				zlog << "Field-Cool..." << std::endl;
			}
			program::field_cool();
			break;

		case 6:
			if(vmpi::my_rank==0){
				std::cout << "Temperature-Pulse..." << std::endl;
				zlog << "Temperature-Pulse..." << std::endl;
			}
			program::temperature_pulse();
			break;
			
		case 7:
			if(vmpi::my_rank==0){
				std::cout << "HAMR-Simulation..." << std::endl;
				zlog << "HAMR-Simulation..." << std::endl;
			}
			program::hamr();
			break;
			
		case 8:
			if(vmpi::my_rank==0){
				std::cout << "CMC-Anisotropy..." << std::endl;
				zlog << "CMC-Anisotropy..." << std::endl;
			}
			program::cmc_anisotropy();
			break;
			
		case 9:
			if(vmpi::my_rank==0){
				std::cout << "Hybrid-CMC..." << std::endl;
				zlog << "Hybrid-CMC..." << std::endl;
			}
			program::hybrid_cmc();
			break;
			
      case 10:
         if(vmpi::my_rank==0){
            std::cout << "Reverse-Hybrid-CMC..." << std::endl;
            zlog << "Reverse-Hybrid-CMC..." << std::endl;
         }
         program::reverse_hybrid_cmc();
         break;

      case 11:
         if(vmpi::my_rank==0){
            std::cout << "LaGrange-Multiplier..." << std::endl;
            zlog << "LaGrange-Multiplier..." << std::endl;
         }
         program::lagrange_multiplier();
         break;

      case 12:
         if(vmpi::my_rank==0){
            std::cout << "partial-hysteresis-loop..." << std::endl;
            zlog << "partial-hysteresis-loop..." << std::endl;
         }
         program::partial_hysteresis_loop();
         break;

		case 50:
			if(vmpi::my_rank==0){
				std::cout << "Diagnostic-Boltzmann..." << std::endl;
				zlog << "Diagnostic-Boltzmann..." << std::endl;
			}
			program::boltzmann_dist();
			break;
		
		default:{
			std::cerr << "Unknown Internal Program ID "<< sim::program << " requested, exiting" << std::endl;
			zlog << "Unknown Internal Program ID "<< sim::program << " requested, exiting" << std::endl;
			exit (EXIT_FAILURE);
			}
	}

	// output Monte Carlo Statistics if applicable
	if(sim::integrator==3 || sim::integrator==4){
		std::cout << "Constrained Monte Carlo Statistics:" << std::endl;
		std::cout << "\tTotal moves: " << cmc::mc_total << std::endl;
		std::cout << "\t" << (cmc::mc_success/cmc::mc_total)*100.0 << "% Accepted" << std::endl;
		std::cout << "\t" << (cmc::energy_reject/cmc::mc_total)*100.0 << "% Rejected (Energy)" << std::endl;
		std::cout << "\t" << (cmc::sphere_reject/cmc::mc_total)*100.0 << "% Rejected (Sphere)" << std::endl;
		
	}

	//program::LLB_Boltzmann();

   // optionally save checkpoint file
   if(sim::save_checkpoint_flag && !sim::save_checkpoint_continuous_flag) save_checkpoint();

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

		case 4: // Hybrid Constrained Monte Carlo
			for(int ti=0;ti<n_steps;ti++){
				sim::ConstrainedMonteCarloMonteCarlo();
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
				terminaltextcolor(RED);
				std::cerr << "Error - Monte Carlo Integrator unavailable for parallel execution" << std::endl;
				terminaltextcolor(WHITE);
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
				terminaltextcolor(RED);
				std::cerr << "Error - Constrained Monte Carlo Integrator unavailable for parallel execution" << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
				// increment time
				increment_time();
			}
			break;
			
		default:{
			terminaltextcolor(RED);
			std::cerr << "Unknown integrator type "<< sim::integrator << " requested, exiting" << std::endl;
			terminaltextcolor(WHITE);
			exit (EXIT_FAILURE);
			}
	}
	
	return EXIT_SUCCESS;
}

} // Namespace sim


