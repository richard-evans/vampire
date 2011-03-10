///
/// @file
/// @brief Contains the Curie Temperature program
///
/// @details Performs a temperature loop to calculate temperature dependent magnetisation
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	09/03/2011
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

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
/// @version 1.1
/// @date    09/03/2011
///
/// @param[in] init Determines whether the system is initialised randomly (0) or ordered (1)
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		11/01/2010
///	Revision:	09/03/2011
///=====================================================================================
///
int curie_temperature(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::curie_temperature has been called" << std::endl;}

	// Set starting temperature
	sim::temperature=sim::Tmin;

	// Perform Temperature Loop
	while(sim::temperature<=sim::Tmax){
				
		// (re)-initialise mean m and counter (mean ought to be a vector[mat]...)
		double meanT=0.0;
		double count=0.0;

		// Equilibrate system
		sim::integrate(sim::equilibration_time);
		
		// Simulate system
		while(sim::time<sim::loop_time){
			
			// Integrate system
			sim::integrate(sim::partial_time);
		
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
		
		// Increment temperature
		sim::temperature+=sim::delta_temperature;
		
	} // End of temperature loop


		
	return EXIT_SUCCESS;
}

}//end of namespace program
