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
/// @brief Contains the LLG (Heun) integrator
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

//Function prototypes
int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);

namespace LLG_arrays{

	// Local arrays for LLG integration
	std::vector <double> x_euler_array;
	std::vector <double> y_euler_array;
	std::vector <double> z_euler_array;

	std::vector <double> x_heun_array;
	std::vector <double> y_heun_array;
	std::vector <double> z_heun_array;

	std::vector <double> x_spin_storage_array;
	std::vector <double> y_spin_storage_array;
	std::vector <double> z_spin_storage_array;

	std::vector <double> x_initial_spin_array;
	std::vector <double> y_initial_spin_array;
	std::vector <double> z_initial_spin_array;

	bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}

namespace sim{

/// @brief LLG Initialisation function
///
/// @details Resizes arrays used for Heun integration
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
int LLGinit(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}

	using namespace LLG_arrays;

	x_spin_storage_array.resize(atoms::num_atoms,0.0);
	y_spin_storage_array.resize(atoms::num_atoms,0.0);
	z_spin_storage_array.resize(atoms::num_atoms,0.0);

	x_initial_spin_array.resize(atoms::num_atoms,0.0);
	y_initial_spin_array.resize(atoms::num_atoms,0.0);
	z_initial_spin_array.resize(atoms::num_atoms,0.0);

	x_euler_array.resize(atoms::num_atoms,0.0);
	y_euler_array.resize(atoms::num_atoms,0.0);
	z_euler_array.resize(atoms::num_atoms,0.0);

	x_heun_array.resize(atoms::num_atoms,0.0);
	y_heun_array.resize(atoms::num_atoms,0.0);
	z_heun_array.resize(atoms::num_atoms,0.0);

	LLG_set=true;

  	return EXIT_SUCCESS;
}

/// @brief LLG Heun Integrator Corrector
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
/// @date    07/02/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
int LLG_Heun(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_Heun has been called" << std::endl;}

	using namespace LLG_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGinit();

	// Local variables for system integration
	const int num_atoms = atoms::num_atoms;

	// Store initial spin positions
   for(int atom = 0; atom < num_atoms; ++atom) x_initial_spin_array[atom] = atoms::x_spin_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) y_initial_spin_array[atom] = atoms::y_spin_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) z_initial_spin_array[atom] = atoms::z_spin_array[atom];


	// Calculate fields
	calculate_spin_fields(0,num_atoms);
	calculate_external_fields(0,num_atoms);

	// Calculate Euler Step
	for(int atom=0;atom<num_atoms;atom++){

		const int imaterial=atoms::type_array[atom];
		const double one_oneplusalpha_sq = mp::material[imaterial].one_oneplusalpha_sq; // material specific alpha and gamma
		const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;

		// Store local spin in Sand local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		// Calculate Delta S
		x_euler_array[atom]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
		y_euler_array[atom]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
		z_euler_array[atom]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

   }

   // cast namespaced variable to local constant
   const double mpdt = mp::dt;

   // Calculate Euler Step
   for(int atom = 0; atom < num_atoms; ++atom) x_spin_storage_array[atom] = x_initial_spin_array[atom] + x_euler_array[atom]*mpdt;
   for(int atom = 0; atom < num_atoms; ++atom) y_spin_storage_array[atom] = y_initial_spin_array[atom] + y_euler_array[atom]*mpdt;
   for(int atom = 0; atom < num_atoms; ++atom) z_spin_storage_array[atom] = z_initial_spin_array[atom] + z_euler_array[atom]*mpdt;

   // calculate normalised spin length (holding in heun array_x for vectorizability)
   for(int atom = 0; atom < num_atoms; ++atom){
      const double x = x_spin_storage_array[atom];
      const double y = y_spin_storage_array[atom];
      const double z = z_spin_storage_array[atom];
      x_heun_array[atom] = 1.0/sqrt(x*x + y*y + z*z);
   }

   // normalise spin length
   for(int atom = 0; atom < num_atoms; ++atom) atoms::x_spin_array[atom] = x_spin_storage_array[atom] * x_heun_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) atoms::y_spin_array[atom] = y_spin_storage_array[atom] * x_heun_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) atoms::z_spin_array[atom] = z_spin_storage_array[atom] * x_heun_array[atom];

	// Recalculate spin dependent fields
	calculate_spin_fields(0,num_atoms);

	// Calculate Heun Gradients
	for(int atom=0;atom<num_atoms;atom++){

		const int imaterial=atoms::type_array[atom];;
		const double one_oneplusalpha_sq = mp::material[imaterial].one_oneplusalpha_sq;
		const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;

		// Store local spin in Sand local field in H
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
									atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
									atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};

		// Calculate Delta S and store dS in heun array
		x_heun_array[atom]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
		y_heun_array[atom]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
		z_heun_array[atom]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

	}

	// Calculate Heun Step
   const double halfdt = mp::half_dt;
	for(int atom=0;atom<num_atoms;atom++) x_spin_storage_array[atom] = x_initial_spin_array[atom] + halfdt* (x_euler_array[atom]+x_heun_array[atom]);
   for(int atom=0;atom<num_atoms;atom++) y_spin_storage_array[atom] = y_initial_spin_array[atom] + halfdt* (y_euler_array[atom]+y_heun_array[atom]);
   for(int atom=0;atom<num_atoms;atom++) z_spin_storage_array[atom] = z_initial_spin_array[atom] + halfdt* (z_euler_array[atom]+z_heun_array[atom]);

   // calculate normalised spin length (holding in heun array_x for vectorizability)
   for(int atom = 0; atom < num_atoms; ++atom){
      const double x = x_spin_storage_array[atom];
      const double y = y_spin_storage_array[atom];
      const double z = z_spin_storage_array[atom];
      x_heun_array[atom] = 1.0/sqrt(x*x + y*y + z*z);
   }

   // normalise spin length
   for(int atom = 0; atom < num_atoms; ++atom) atoms::x_spin_array[atom] = x_spin_storage_array[atom] * x_heun_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) atoms::y_spin_array[atom] = y_spin_storage_array[atom] * x_heun_array[atom];
   for(int atom = 0; atom < num_atoms; ++atom) atoms::z_spin_array[atom] = z_spin_storage_array[atom] * x_heun_array[atom];

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
int LLG_Heun_cuda(){


	return EXIT_SUCCESS;
}

}
