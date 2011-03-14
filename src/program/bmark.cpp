///
/// @file
/// @brief Contains the standard benchmark program
///
/// @details Simulates a system for a number of timesteps
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
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace program{
	
int bmark(){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::bmark has been called" << std::endl;}

	// Simulate system
	while(sim::time<sim::total_time){
		sim::integrate(sim::partial_time);

		// Calculate mag_m, mag after sim::partial_time steps
      stats::mag_m();

		vout::data();

	} // end of time loop
	
	return EXIT_SUCCESS;
}

}//end of namespace program

