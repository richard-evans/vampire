///
/// @file
/// @brief Contains the Two Temperature Pulse program
///
/// @details Performs a time series with variable temperature
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    06/04/2011
/// @internal
///	Created:		06/04/2011
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
/// @brief Function to calculate a time series with two temperature model heating/cooling profile
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    06/04/2011
/// 
/// @internal
///	Created:		10/03/2010
///	Revision:	06/04/2011
///=====================================================================================
///
void two_temperature_pulse(){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::two_temperature_pulse has been called" << std::endl;}
	
	const double Ce = 7.0E02; //electron specific heat 
	const double Cl = 3.0E06; //phonon specific heat
	const double G = 17.0E17 ;//electron coupling constant
	const double Q = 1.0E10; // heatsink specific heat
	const double K = 20.0E17; // phonon-heatsink coupling
	
	//double pump_time=20.0e-15; // Seconds //moved to sim
	//double pump_power=2.4e22; // ? // moved to sim
	
	// Set equilibration temperature and field
	sim::temperature=sim::Teq;

	// Equilibrate system
	while(sim::time<sim::equilibration_time){
		
		sim::integrate(sim::partial_time);
		
		// Calculate magnetisation statistics
		stats::mag_m();
		
		// Output data
		vout::data();
	}
	
	double Te=sim::temperature;
	double Tp=sim::temperature;
	
	//loop 5 times
	for(int r=0; r<sim::runs;r++){
		
	// record starting time after equiibration/last pulse
	int start_time=sim::time;

	// Simulate temperature pulse
	while(sim::time<sim::total_time+start_time){

		// Calculate Temperature
		double time_from_start=mp::dt_SI*double(sim::time-start_time);
		double pump=sim::pump_power*exp(-((time_from_start-3.*sim::pump_time)/(sim::pump_time) )*((time_from_start-3.*sim::pump_time)/(sim::pump_time) ));

		Te = (-G*(Te-Tp)+pump)*mp::dt_SI/(Ce*Te) + Te;
		Tp = ( G*(Te-Tp)     )*mp::dt_SI/Cl + Tp; // - (Tp-sim::Tmin)*3.0e11*mp::dt_SI;
		
		sim::temperature=Te;
		vinfo << sim::time << "\t" << Te << "\t" << Tp << std::endl;
		
		// Integrate system
		sim::integrate(sim::partial_time);
		
		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();

	}
	
	}

} // end of two temperature pulse

} // end of namespace program

