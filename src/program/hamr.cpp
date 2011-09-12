///
/// @file
/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@nanohpc.com
/// @version 1.0
/// @date    10/06/2011
/// @internal
///	Created:		10/06/2011
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

/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
/// @details Performs a time series with moving localised head field and temperature pulse
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    10/06/2011
///
/// @internal
///	Created:		10/06/2011
///	Revision:	--
///=====================================================================================
///
void hamr(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hamr has been called" << std::endl;}

		// Set equilibration temperature and field
		sim::temperature=sim::Teq;
		
		// Disable laser
		sim::head_laser_on=false;

		// Equilibrate system
		while(sim::time<sim::equilibration_time){
			
			sim::integrate(sim::partial_time);
			
			// Calculate magnetisation statistics
			stats::mag_m();
			
			// Output data
			vout::data();
		}

		// now enable laser
		sim::head_position[0]=0.0;
		sim::head_position[1]=cs::system_dimensions[1]*0.5; // A
		sim::head_speed=30.0; // nm/ns
		sim::head_laser_on=true;
		
		int start_time=sim::time;
		
		// Perform HAMR
		while(sim::time<sim::total_time+start_time){
			
			// Integrate system
			sim::integrate(sim::partial_time);
			
			// Calculate magnetisation statistics
			stats::mag_m();

			// Output data
			vout::data();

		}
	
} // end of hamr()

}//end of namespace program
