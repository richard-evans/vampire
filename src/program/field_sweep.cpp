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
#include <sstream>
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
int field_sweep(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hysteresis has been called" << std::endl;}


	//reads in a file called field-sweep
	//this fiel should be of the format:
	// Hstart , Hend, increment, n_steps
	//this can have as many lines as you like
	//they will be run in the order of the file


	std::ifstream parameter_file;
  parameter_file.open("field-sweep");

	if (!parameter_file) {
	    std::cout << "Unable to open file field-sweep";
			err::vexit();   // call system to stop
	}


	double Hstart, Hend;
	int n_steps;
	double increment;
	int start_time;
  std::string line; // declare a string to hold line of text
  while(getline(parameter_file,line)){
    std::stringstream line_stream(line);
		if (line.length() != 0 && line[0] != '#'){
	//	std::cout << line << std::endl;
    line_stream >> Hstart >> Hend >> increment >> n_steps;
	//	std::cout << Hstart << '\t' << Hend << '\t' << increment << '\t' << n_steps <<std::endl;

		if (increment < 0) increment = -increment;
		if (Hstart > Hend) increment = -increment;


	// Equilibrate system in saturation field, i.e. the largest between equilibration and maximum field set by the user
   if(sim::Heq >= sim::Hmax){
	   sim::H_applied=sim::Heq;
   }
   else{
   	sim::H_applied=sim::Hmax;
   }

	// Initialise sim::integrate only if it not a checkpoint
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else sim::integrate(sim::equilibration_time);

   // Hinc must be positive
	//int iHinc=vmath::iround(double(fabs(sim::Hinc))*1.0E6);
	double field = Hstart;
	bool run = true;
		// Perform Field loop-field
		while (run){

			// Set applied field (Tesla)
			sim::H_applied=field;

			// Reset start time
			start_time=sim::time;

			// Reset mean magnetisation counters
			stats::reset();

			sim::loop_time = n_steps;
			// Integrate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);

				// Calculate mag_m, mag
				stats::update();

			} // End of integration loop
//			sim::iH=int64_t(Hfield); //sim::iH+=iHinc;

			// Output to screen and file after each field
			vout::data();
			if (Hstart > Hend && field <= Hend) run = false;
			if (Hstart < Hend && field >= Hend) run = false;
			field = field + increment;


		} // End of field loop

		}
	}
	parameter_file.close();
	return EXIT_SUCCESS;

} // End of hysteresis program

}//end of namespace program
