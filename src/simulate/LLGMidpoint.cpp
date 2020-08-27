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
/// @brief Contains the LLG (Midpoint) integrator
///
/// @details The LLG integrator...
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
/// @date    07/02/2011
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
#include "LLG.hpp"
#include "material.hpp"
#include "sim.hpp"


namespace sim{

/// @brief LLG Midpoint Integrator
///
/// @callgraph
/// @callergraph
///
/// @details Integrates the system using the LLG and Heun solver
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    14/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int LLG_Midpoint(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_Midpoint has been called" << std::endl;}

	using namespace LLG_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGinit();

	// Local variables for system integration
	const int num_atoms=atoms::num_atoms;

	// Store initial spin positions
	for(int atom=0;atom<num_atoms;atom++){
		x_initial_spin_array[atom] = atoms::x_spin_array[atom];
		y_initial_spin_array[atom] = atoms::y_spin_array[atom];
		z_initial_spin_array[atom] = atoms::z_spin_array[atom];
	}

	// Calculate fields
	calculate_spin_fields(0,num_atoms);
	calculate_external_fields(0,num_atoms);

	// Calculate Predictor Step
	for(int atom=0;atom<num_atoms;atom++){

		const int imaterial=atoms::type_array[atom];
		const double alpha = mp::material[imaterial].alpha;
		const double beta  = -1.0*mp::dt*mp::material[imaterial].one_oneplusalpha_sq*0.5;
		const double beta2 = beta*beta;

		// Store local spin in S and local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		// Calculate F = [H + alpha* (S x H)]
		const double F[3] = {H[0] + alpha*(S[1]*H[2]-S[2]*H[1]),
									H[1] + alpha*(S[2]*H[0]-S[0]*H[2]),
									H[2] + alpha*(S[0]*H[1]-S[1]*H[0])};

		const double FdotF = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
		const double beta2FdotS = beta2*(F[0]*S[0] + F[1]*S[1] + F[2]*S[2]);
		const double one_o_one_plus_beta2FdotF = 1.0/(1.0 + beta2*FdotF);
		const double one_minus_beta2FdotF = 1.0 - beta2*FdotF;

		// Calculate intermediate spin position (S + S')/2
		x_spin_storage_array[atom] = (S[0] + one_o_one_plus_beta2FdotF*(S[0]*one_minus_beta2FdotF + 2.0*(beta*(F[1]*S[2]-F[2]*S[1]) + F[0]*beta2FdotS)))*0.5;
		y_spin_storage_array[atom] = (S[1] + one_o_one_plus_beta2FdotF*(S[1]*one_minus_beta2FdotF + 2.0*(beta*(F[2]*S[0]-F[0]*S[2]) + F[1]*beta2FdotS)))*0.5;
		z_spin_storage_array[atom] = (S[2] + one_o_one_plus_beta2FdotF*(S[2]*one_minus_beta2FdotF + 2.0*(beta*(F[0]*S[1]-F[1]*S[0]) + F[2]*beta2FdotS)))*0.5;

	}

	// Store Midpoint to spin array
	for(int atom=0;atom<num_atoms;atom++){
		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
		atoms::z_spin_array[atom]=z_spin_storage_array[atom];
	}

	// Recalculate spin dependent fields
	calculate_spin_fields(0,num_atoms);

	// Calculate Corrector Step
	for(int atom=0;atom<num_atoms;atom++){

		const int imaterial=atoms::type_array[atom];
		const double alpha = mp::material[imaterial].alpha;
		const double beta  = -1.0*mp::dt*mp::material[imaterial].one_oneplusalpha_sq*0.5;
		const double beta2 = beta*beta;

		// Store local spin in S and local field in H
		const double M[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double S[3] = {x_initial_spin_array[atom],y_initial_spin_array[atom],z_initial_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		// Calculate F = [H + alpha* (M x H)]
		const double F[3] = {H[0] + alpha*(M[1]*H[2]-M[2]*H[1]),
									H[1] + alpha*(M[2]*H[0]-M[0]*H[2]),
									H[2] + alpha*(M[0]*H[1]-M[1]*H[0])};

		const double FdotF = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
		const double beta2FdotS = beta2*(F[0]*S[0] + F[1]*S[1] + F[2]*S[2]);
		const double one_o_one_plus_beta2FdotF = 1.0/(1.0 + beta2*FdotF);
		const double one_minus_beta2FdotF = 1.0 - beta2*FdotF;

		// Calculate final spin position
		atoms::x_spin_array[atom] = one_o_one_plus_beta2FdotF*(S[0]*one_minus_beta2FdotF + 2.0*(beta*(F[1]*S[2]-F[2]*S[1]) + F[0]*beta2FdotS));
		atoms::y_spin_array[atom] = one_o_one_plus_beta2FdotF*(S[1]*one_minus_beta2FdotF + 2.0*(beta*(F[2]*S[0]-F[0]*S[2]) + F[1]*beta2FdotS));
		atoms::z_spin_array[atom] = one_o_one_plus_beta2FdotF*(S[2]*one_minus_beta2FdotF + 2.0*(beta*(F[0]*S[1]-F[1]*S[0]) + F[2]*beta2FdotS));
	}

	return EXIT_SUCCESS;
}

/// @brief LLG Heun Integrator (CUDA)
///
/// @callgraph
/// @callergraph
///
/// @details Integrates the system using the LLG and Heun solver
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
int LLG_Midpoint_cuda(){


	return EXIT_SUCCESS;
}

}
