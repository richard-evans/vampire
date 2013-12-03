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
/// @brief Contains the Hysteresis program
///
/// @details Performs a field loop to calculate field dependent magnetisation
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    21/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	21/03/2011
///=====================================================================================
///

// Standard Libraries
#include <cstdlib>

// Vampire Header files
#include "vmath.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"


namespace program{

/// @brief Function to calculate the hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// ordered. After initialisation a whole hysteresis loop of the system and coercivity are calculated.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Weijia Fan, wf507@york.ac.uk
/// @version 1.0
/// @date    27/01/2010
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		27/01/2010
///	Revision:	  ---
///=====================================================================================
///
int hysteresis(){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hysteresis has been called" << std::endl;}
	
	// Equilibrate system in saturation field
	sim::H_applied=sim::Heq;
	sim::integrate(sim::equilibration_time);
		
	// Setup min and max fields and increment (uT)
	int iHmax=vmath::iround(double(sim::Hmax)*1.0E6);
	int iHmin=vmath::iround(double(sim::Hmin)*1.0E6);
   // Hinc must be positive
	int iHinc=vmath::iround(double(fabs(sim::Hinc))*1.0E6);

	// Perform Field Loop
	for(int parity=-1;parity<2;parity+=2){
		for(int H=-iHmax;H<=iHmax;H+=iHinc){
			
			// Set applied field (Tesla)
			sim::H_applied=double(H)*double(parity)*1.0e-6;
			
			// Reset start time
			int start_time=sim::time;
			
			// Reset mean magnetisation counters
			stats::mag_m_reset();

			// Integrate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);
			
				// Calculate mag_m, mag
				stats::mag_m();

			}

			// Output to screen and file after each field
			vout::data();
			
		} // End of field loop
	} // End of parity loop

	return EXIT_SUCCESS;
  }

}//end of namespace program


