#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include <iostream>

namespace program{

/// @brief Function to calculate the temperature dependence of the magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// accoring to the input flag - either randomly or ordered.For the ordered case the temperature sequence
/// increases from zero, for the random case the temperature decreases from the maximum temperature. After
/// initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
/// @section notes Implementation Notes
/// Capable of hot>cold or cold>hot calculation. 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    11/01/2010
///
/// @param[in] init Determines whether the system is initialised randomly (0) or ordered (1)
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
int curie_temperature(bool init){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::curie_temperature has been called" << std::endl;}
	
	// function prototypes

	//int output_pov_file();

	// Setup LLG arrays
	sim::initialise();

	// Initialise spins to random state
	for(int atom =0;atom<atoms::num_atoms;atom++){
		if(init==false){
			atoms::x_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			atoms::y_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			atoms::z_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
		}
		else{
		    double parity = 1.0; //double(atoms::grain_array[atom]%2);
			atoms::x_spin_array[atom]=0.01;			
			atoms::y_spin_array[atom]=0.0;			
			atoms::z_spin_array[atom]=2.0*parity-1.0;
		}
	}
		//sim::H_applied=0.0;
	// Set up loop variables
	
	sim::constraint_phi=0.0;
	sim::constraint_theta=0.0;
	
	sim::CMCinit();

	//vout::pov_file();
	sim::integrator=3;

	//      Perform Temperature Loop
	for(int temperature=0;temperature<=1500;temperature+=10){
		// Set system temperature
		sim::temperature=double(temperature); 
		
		double meanT=0.0;
		double count=0.0;

		// Equilibrate system
		sim::integrate(sim::equilibration_time);
		//sim::LLG(sim::equilibration_time);
		
		// Simulate system
		for(sim::time=0;sim::time<sim::loop_time;sim::time+=sim::partial_time){
			
			// Integrate system
			sim::integrate(sim::partial_time);
		
			// Calculate LLG
			//sim::LLG(sim::partial_time);
			
			// Calculate mag_m, mag
			stats::mag_m();
			meanT+=stats::total_mag_m_norm;
			count+=1.0;
		}
		// Output to screen and file after each temperature
		if(vmpi::my_rank==0){
			std::cout << sim::temperature << "\t" << stats::total_mag_m_norm;
			std::cout << "\t" << stats::total_mag_norm[0];
			std::cout << "\t" << stats::total_mag_norm[1];
			std::cout << "\t" << stats::total_mag_norm[2];
			std::cout << "\t" << meanT/count;
			std::cout << std::endl;
			vmag << sim::temperature << "\t" << stats::total_mag_m_norm << "\t" << meanT/count << std::endl;
		}
		
		//vout::pov_file();
		
	} // End of temperature loop


		
	return EXIT_SUCCESS;
}

}//end of namespace program
