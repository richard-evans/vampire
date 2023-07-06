//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
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
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void time_series(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::time_series has been called" << std::endl;

	double temp=sim::temperature;

   // Set equilibration temperature only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
	   // Set equilibration temperature
	   sim::temperature=sim::Teq;
   }

	// Output data
	vout::data();

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();
	}

   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set simulation temperature
	   sim::temperature = temp;

      // Reset mean magnetisation counters
      stats::reset();

   }

	// Perform Time Series
	while(sim::time<sim::equilibration_time+sim::total_time){

		// Integrate system
		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

	}

}

}//end of namespace program
