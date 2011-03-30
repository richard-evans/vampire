///
/// @file
/// @brief Contains the Field Cooling program
///
/// @details Performs a time series with decreasing temperature
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Function to calculate the a field cooled magnetisation
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
///
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///
void field_cool(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::field_cool has been called" << std::endl;}

	// Set equilibration temperature and field
	sim::temperature=sim::Teq;
	
	// Equilibrate system
	while(sim::time<sim::equilibration_time){
		
		sim::integrate(sim::partial_time);
		
		// Output data
		vout::data();
	}
	
	int start_time=sim::time;
	
	// Perform Field Cooling
	while(sim::time<sim::total_time+start_time){

		// Calculate Temperature
		double time_from_start=mp::dt_SI*double(sim::time-start_time);
		if(sim::cooling_function_flag==0){
			sim::temperature = sim::Tmin + (sim::Tmax-sim::Tmin)*exp(-time_from_start/sim::cooling_time);
		}
		else{
			sim::temperature = sim::Tmin + (sim::Tmax-sim::Tmin)*exp(-(time_from_start)*(time_from_start)/((sim::cooling_time)*(sim::cooling_time)));
		}
		
		// Integrate system
		sim::integrate(sim::partial_time);
		
		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();

	}

}

}//end of namespace program
