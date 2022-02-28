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
#include "montecarlo.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Function to calculate the temperature dependence of the anisotropy and magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature, where the constraint angles
/// are cycled. The system is initialised all spins along the constraint direction. After initialisation
/// the sytem is equilibrated for sim::equilibration timesteps before statistics are collected.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    15/09/2011
///
/// @return void
///
/// @internal
///	Created:		15/09/2011
///	Revision:	--/--/----
///=====================================================================================
///
void cmc_anisotropy(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::cmc_anisotropy has been called" << std::endl;}

	// Check integrator is CMC, if not then exit disgracefully
	if(sim::integrator!=3){
		err::zexit("Program CMC-anisotropy requires Constrained Monte Carlo as the integrator. Check input file.");
	}

	// set minimum rotational angle
   // if checkpoint is loaded, then update minimum values of temperature and constrained angles
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
      zlog << zTs() << "Constrained angles theta loaded from checkpoint " << sim::constraint_theta <<  std::endl;
   }
   else sim::constraint_theta=sim::constraint_theta_min;

	// perform rotational angle sweep
	while(sim::constraint_theta<=sim::constraint_theta_max){

		// set minimum azimuthal angle
      // if checkpoint is loaded, then update minimum values of temperature and constrained angles
	   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
         zlog << zTs() << "Constrained angle phi loaded from checkpoint: " << sim::constraint_phi << std::endl;
      }
		else sim::constraint_phi=sim::constraint_phi_min;

		// perform azimuthal angle sweep
		while(sim::constraint_phi<=sim::constraint_phi_max){

			// Re-initialise spin moments for CMC
			montecarlo::CMCinit();

			// Set starting temperature
         // if checkpoint is loaded, then update minimum values of temperature
	      if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
            //zlog << zTs() << "Constrained temperature loaded from checkpoint: " << sim::temperature << "K and flag = " << sim::checkpoint_loaded_flag <<std::endl;
            // In order to correctly restart the simulation, the flag that chech if the checkpoint
            // load is set to restart or load needs to be silenced
            sim::checkpoint_loaded_flag=false;
            sim::temperature += sim::delta_temperature;
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

			// Increment azimuthal angle
			sim::constraint_phi+=sim::constraint_phi_delta;
			sim::constraint_phi_changed=true;

		} // End of azimuthal angle sweep
      if(vout::gnuplot_array_format) zmag << std::endl;

		// Increment rotational angle
		sim::constraint_theta+=sim::constraint_theta_delta;
		sim::constraint_theta_changed=true;

	} // End of rotational angle sweep

	return;
}

}//end of namespace program
