//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2025. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

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

//------------------------------------------------------------------------------
// Function to cool the system under the action of a constant applied field at
// a defined cooling rate in K/ns
//------------------------------------------------------------------------------
void hamr_cool(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hamr_cool has been called" << std::endl;}



	// Set equilibration temperature and field
	sim::temperature = sim::Teq;

	//--------------------------------------------------------------
	// Equilibrate system
	//--------------------------------------------------------------
	while( sim::time < sim::equilibration_time ){

		// Integrate system
		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

	}

	// Cooling and field rate parameters
	const double cooling_rate = 200.0; // K/ns
	const double field_rate = 8.5333; // T/ns

	// derived delta T and delta B values per timestep
	const double dT_per_dt = cooling_rate * 1.0e9 * mp::dt_SI; // K per timestep
	const double dB_per_dt = field_rate   * 1.0e9 * mp::dt_SI;   // T per timestep
	const double dT_max = sim::Tmax - sim::Tmin;

	const double B_max = sim::H_applied;
	const double dB_max = -2.0 * B_max;

	// Calculate number of timesteps to reach final temperature and field
	const int cooling_time  = dT_max / dT_per_dt;
	const int field_sw_time = dB_max / dB_per_dt;

	std::cout << dT_per_dt << "\t" << dB_per_dt << "\t" << dT_max << "\t" << dB_max << std::endl;

	// Print information about simulation time to user
	std::cout << "Time steps to reach final temperature: " << cooling_time << std::endl;
	if( cooling_time > sim::total_time){
		std::cerr << "WARNING: specified time steps of " << sim::total_time << " is insufficient to reach final temperature during hamr cooling program" << std::endl;
		zlog << zTs() << "WARNING: specified time steps of " << sim::total_time << " is insufficient to reach final temperature during hamr cooling program" << std::endl;
	}

	// save starting time after equilibration
	uint64_t start_time = sim::time;

	//--------------------------------------------------------------
	// Perform field cooling
	//--------------------------------------------------------------
	while( sim::time < sim::total_time + start_time ){

      // loop over partial time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate dynamic temperature and field
         double time_from_start = double( sim::time - start_time );

         double Tlinear = sim::Tmax - (dT_max)*(time_from_start/cooling_time);
         if( Tlinear >= sim::Tmin ) sim::temperature = Tlinear;
         else sim::temperature = sim::Tmin;

			double B_linear = B_max - (dB_max)*(time_from_start/field_sw_time);
         if( B_linear >= -B_max ) sim::H_applied = B_linear;
         else sim::H_applied = -B_max;

         // Integrate system
         sim::integrate(1);

		}

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

		}

} // end of hamr_cool()

} //end of namespace program
