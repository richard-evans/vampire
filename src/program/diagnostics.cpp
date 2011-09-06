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
#include "vmath.hpp"


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

void boltzmann_dist(){
	
	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::boltzmann_dist has been called" << std::endl;

	// array for binning spin angle
	std::vector<double> bin(181,0.0);
	
	// Set starting temperature
	sim::temperature=sim::temperature;

	// Equilibrate system
	sim::integrate(sim::equilibration_time);
	
	// Simulate system
	while(sim::time<sim::total_time+sim::equilibration_time){
		
		// Integrate system
		sim::integrate(sim::partial_time);
		
		// Calculate magnetisation statistics
		for(int atom=0; atom<atoms::num_atoms; atom++){
			double angle = acos(atoms::z_spin_array[atom])*180.0/M_PI;
			double id = vmath::iround(angle+0.5);
			bin[id]+=1.0;
		}
	}
	
	// Find max probability and max P
	double maxp=0.0;
	double maxP=0.0;
	for(int b=0;b<181;b++){
		double energy = mp::material[0].Ku1_SI;
		double P = sin(double (b)*M_PI/180)*exp((energy*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0))/(sim::temperature*1.3806503e-23));
		if((bin[b])>maxp) maxp=bin[b];
		if(P>maxP) maxP=P;
	}

	// Output data
	vmag << "# " << mp::material[0].Ku1_SI/(sim::temperature*1.3806503e-23) << std::endl;
	for(int b=0;b<181;b++){
		double energy = mp::material[0].Ku1_SI;
		double P = sin(double (b)*M_PI/180)*exp((energy*sin(double (b)*M_PI/180.0)*sin(double (b)*M_PI/180.0))/(sim::temperature*1.3806503e-23));
		vmag << b << "\t" << (bin[b]+bin[180-b])/(2.0*maxp) << "\t" << P/maxP << std::endl;
	}
	
}

}//end of namespace program

