///
/// @file
/// @brief Contains the Monte Carlo integrator
///
/// @details The Monte Carlo integrator...
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation. 
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"

namespace sim{
	
/// @brief Monte Carlo Integrator
///
/// @callgraph
/// @callergraph
///
/// @details Integrates the system using a Monte Carlo solver with tuned step width 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int MonteCarlo(int istart, int iend){
	
	// Check for calling of function
	if(err::check==true) std::cout << "sim::MonteCarlo has been called" << std::endl;

	// calculate number of seps to calculate
	int nmoves = iend-istart;
	if(nmoves <=0){
		std::cerr << "Error in Monte Carlo Integrator, start and end ranges must satisfy the following:" << std::endl;
		std::cerr << "      istart  <  iend" << std::endl;
		std::cerr << "      istart  >= 0" << std::endl;
		std::cerr << "      iend    >  0" << std::endl;
		std::cerr << "      istart is currently: " << istart << std::endl;
		std::cerr << "      iend is currently:   " << iend << std::endl;
		std::cerr << "Exiting" << std::endl;
		//vexit(vmpi::my_rank);
	}

// Declare arrays for spin statesSave spin position
		double Sold[3];
		double Snew[3];
		
	// Calculate range for move
	double sigma=1.0;
	
	// Temporaries
	int atom=0;
	double r=1.0;
	double Eold=0.0;
	double Enew=0.0;
	double DE=0.0;
	double kBTBohr = sim::temperature*1.3806503e-23/9.27400915e-24;
	
	// loop over natoms to form a single Monte Carlo step
	for(int i=0;i<nmoves; i++){
		
		// pick atom
		atom = istart+int(nmoves*mtrandom::grnd());

		// Save old spin position
		Sold[0] = atoms::x_spin_array[atom];
		Sold[1] = atoms::y_spin_array[atom];
		Sold[2] = atoms::z_spin_array[atom];

		// Calculate new spin position cf Pierre Asselin
		Snew[0] = mtrandom::gaussian()*sigma+Sold[0];
		Snew[1] = mtrandom::gaussian()*sigma+Sold[1];
		Snew[2] = mtrandom::gaussian()*sigma+Sold[2];

		// Calculate new spin length and normalise
		r = 1.0/sqrt (Snew[0]*Snew[0]+Snew[1]*Snew[1]+Snew[2]*Snew[2]); 

		Snew[0]*=r;
		Snew[1]*=r;
		Snew[2]*=r;

		// Calculate current energy
		//Eold = calculate_spin_energy(atom);
		
		// Copy new spin position
		atoms::x_spin_array[atom] = Snew[0];
		atoms::y_spin_array[atom] = Snew[1];
		atoms::z_spin_array[atom] = Snew[2];

		// Calculate new energy
		//Enew = calculate_spin_energy(atom);
		
		// Calculate difference
		DE = Enew-Eold;
		
		// Check for lower energy state and accept unconditionally
		if(DE<0) continue;
		// Otherwise evaluate probability for move
		else{
			if(exp(-DE/kBTBohr) >= mtrandom::grnd()) continue;
			// If rejected reset spin coordinates and continue
			else{
				atoms::x_spin_array[atom] = Sold[0];
				atoms::y_spin_array[atom] = Sold[1];
				atoms::z_spin_array[atom] = Sold[2];
				continue;
			}
		}
	}
	
	return EXIT_SUCCESS;
}

} // End of namespace sim

