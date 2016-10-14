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
#include "material.hpp"
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
int MonteCarlo(){

	// Check for calling of function
	if(err::check==true) std::cout << "sim::MonteCarlo has been called" << std::endl;

	// calculate number of steps to calculate
	int nmoves = atoms::num_atoms;

	// Declare arrays for spin states
	std::valarray<double> Sold(3);
	std::valarray<double> Snew(3);

	// Temporaries
	int atom=0;
	double Eold=0.0;
	double Enew=0.0;
	double DE=0.0;
	const int AtomExchangeType=atoms::exchange_type;

   // Material dependent temperature rescaling
   std::vector<double> rescaled_material_kBTBohr(mp::num_materials);
   std::vector<double> sigma_array(mp::num_materials); // range for tuned gaussian random move
   for(int m=0; m<mp::num_materials; ++m){
      double alpha = mp::material[m].temperature_rescaling_alpha;
      double Tc = mp::material[m].temperature_rescaling_Tc;
      double rescaled_temperature = sim::temperature < Tc ? Tc*pow(sim::temperature/Tc,alpha) : sim::temperature;
      rescaled_material_kBTBohr[m] = 9.27400915e-24/(rescaled_temperature*1.3806503e-23);
      sigma_array[m] = rescaled_temperature < 1.0 ? 0.02 : pow(1.0/rescaled_material_kBTBohr[m],0.2)*0.08;
   }

   double statistics_moves = 0.0;
   double statistics_reject = 0.0;

	// loop over natoms to form a single Monte Carlo step
	for(int i=0;i<nmoves; i++){

      // add one to number of moves counter
      statistics_moves+=1.0;

		// pick atom
		atom = int(nmoves*mtrandom::grnd());

		// get material id
		const int imaterial=atoms::type_array[atom];

      // Calculate range for move
      sim::mc_delta_angle=sigma_array[imaterial];

		// Save old spin position
		Sold[0] = atoms::x_spin_array[atom];
		Sold[1] = atoms::y_spin_array[atom];
		Sold[2] = atoms::z_spin_array[atom];

      // Make Monte Carlo move
      sim::mc_move(Sold, Snew);

		// Calculate current energy
		Eold = sim::calculate_spin_energy(atom, AtomExchangeType);

		// Copy new spin position
		atoms::x_spin_array[atom] = Snew[0];
		atoms::y_spin_array[atom] = Snew[1];
		atoms::z_spin_array[atom] = Snew[2];

		// Calculate new energy
		Enew = sim::calculate_spin_energy(atom, AtomExchangeType);

		// Calculate difference in Joules/mu_B
		DE = (Enew-Eold)*mp::material[imaterial].mu_s_SI*1.07828231e23; //1/9.27400915e-24

		// Check for lower energy state and accept unconditionally
		if(DE<0) continue;
		// Otherwise evaluate probability for move
		else{
			if(exp(-DE*rescaled_material_kBTBohr[imaterial]) >= mtrandom::grnd()) continue;
			// If rejected reset spin coordinates and continue
			else{
				atoms::x_spin_array[atom] = Sold[0];
				atoms::y_spin_array[atom] = Sold[1];
				atoms::z_spin_array[atom] = Sold[2];
            // add one to rejection counter
            statistics_reject += 1.0;
				continue;
			}
		}
	}

   // Save statistics to sim namespace variable
   sim::mc_statistics_moves += statistics_moves;
   sim::mc_statistics_reject += statistics_reject;

	return EXIT_SUCCESS;
}

} // End of namespace sim
