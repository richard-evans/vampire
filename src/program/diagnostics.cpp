///
/// @file
/// @brief Contains the standard benchmark program
///
/// @details Simulates a system for a number of timesteps
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
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"


namespace program{
	
int timestep_scaling(){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::timestep_scaling has been called" << std::endl;}

	std::cout << " Diagnostic - Timestep Scaling " << std::endl;

	// loop over timesteps                                                                                                                                                     
	for(int powerv=18; powerv > 13; powerv--){
		for(int value=1;value<10;value++){
			mp::dt_SI=double(value)*pow(10.0,-1.0*powerv);

			int timesteps = 5.0e-12/mp::dt_SI;

			std::cout << timesteps << std::endl;

			// reset derived parameters                                                                                                                                
			mp::set_derived_parameters();

			double sx = 0.01;
			double sy = 0.0;
			double sz = 1.0;

			double modS = 1.0/sqrt(sx*sx + sy*sy + sz*sz);

			sx*=modS;
			sy*=modS;
			sz*=modS;

			for(int atom=0; atom<atoms::num_atoms; atom++){
				atoms::x_spin_array[atom] = sx;
				atoms::y_spin_array[atom] = sy;
				atoms::z_spin_array[atom] = sz;
			}

			sim::integrate(timesteps);
			stats::mag_m_reset();
			int start_time=sim::time;
			// Simulate system                                                                                                                                         
			while(sim::time<timesteps+start_time){
				sim::integrate(1);
				
				// Calculate mag_m, mag after sim::partial_time steps                                                                                              
				stats::mag_m();

			} // end of time loop                                                                                                                                      
			vmag << mp::dt_SI << "\t";
			std::cout << mp::dt_SI << "\t";
			vout::data();
			}
		}

	return EXIT_SUCCESS;
}

}//end of namespace program

