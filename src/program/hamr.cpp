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
///	Revision:	Andrea Meo, Feb 2022
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "hamr.hpp"
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
///	Revision:	Andrea Meo, 2022
///=====================================================================================
///
	void hamr(){

		// check calling of routine if error checking is activated
		if(err::check==true){std::cout << "program::hamr has been called" << std::endl;}

			// Set partial_time (sim:time-steps-increment) to 1 to ensure good simulation
			sim::partial_time = 1;

			// Set equilibration temperature and field
			sim::temperature=sim::Teq;
		
			// Disable laser
			hamr::head_laser_on=false;

			// Equilibrate system
			while(sim::time<sim::equilibration_time){
			
				sim::integrate(sim::partial_time);
			
				// Calculate magnetisation statistics
				stats::update();
			
				// Output data
				vout::data();
			}

			// now enable laser
			hamr::head_laser_on=true;
			std::cout << " Equilibration terminated, enabling laser ..." << std::endl;
			zlog << zTs() << "Equilibration terminated, enabling laser ..." << std::endl;

			// Set system temperature as minimum temperature
			sim::temperature=sim::Tmin;

			// Reset mean magnetisation counters
			stats::reset();

			// Perform harm continuous simulation
			hamr::hamr_continuous();

			// //-------------------------------------------------------------------------//
			// // force outputting at end of simulation after laser has been switched off
			// //-------------------------------------------------------------------------//
			// // Disable laser
			// hamr::head_laser_on=false;
			// // Switch off external field 
			// sim::H_applied = 0.0;
			// // Set system temperature as minimum temperature
			// sim::temperature=sim::Tmin;
			// std::cout << " Disabling laser and external field and integrating system for 1 time-step" << std::endl;
			// zlog << zTs() << " Disabling laser and external field and integrating system for 1 time-step" << std::endl;
			// // Integrate
			// sim::integrate(sim::partial_time);
			// // Calculate magnetisation statistics
			// stats::update();
			// // Output data
			// vout::data();
			// std::cout << " Outputting system at the end of HAMR simulations\n" << std::endl;
			// zlog << zTs() << "Outputting system at the end of HAMR simulations" << std::endl;

	} // end of hamr()

}//end of namespace program
