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
   // Initialise sim::temperature
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){
      sim::temperature+=sim::delta_temperature;
   }
   else sim::temperature=sim::Tmin;

	// Perform Temperature Loop
	while(sim::temperature<=sim::Tmax){

		// Equilibrate system
		sim::integrate(sim::equilibration_time);

		// Reset mean magnetisation counters
		stats::reset();

		// Reset start time
		int start_time=sim::time;

		// Simulate system
		while(sim::time<sim::loop_time+start_time){

			// Integrate system
			sim::integrate(sim::partial_time);

			// Calculate magnetisation statistics
			stats::update();

		}

		// Output data
		vout::data();

		// Increment temperature
		sim::temperature+=sim::delta_temperature;

	} // End of temperature loop

	return EXIT_SUCCESS;
}

}//end of namespace program
