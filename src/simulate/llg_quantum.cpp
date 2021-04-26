//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

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

namespace LLGQ_arrays{

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
int LLGQinit(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}

	using namespace LLGQ_arrays;

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

namespace internal{

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
void llg_quantum_step(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_Heun has been called" << std::endl;}

	using namespace LLGQ_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGQinit();

	// Local variables for system integration
	const int num_atoms=atoms::num_atoms;
	double xyz[3];		// Local Delta Spin Components
	double S_new[3];	// New Local Spin Moment
	double mod_S;		// magnitude of spin moment

	// Store initial spin positions
	for(int atom=0;atom<num_atoms;atom++){
		x_initial_spin_array[atom] = atoms::x_spin_array[atom];
		y_initial_spin_array[atom] = atoms::y_spin_array[atom];
		z_initial_spin_array[atom] = atoms::z_spin_array[atom];
	}

	// Calculate fields
	sim::calculate_spin_fields(0,num_atoms);
	sim::calculate_external_fields(0,num_atoms);

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
		xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
		xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
		xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

		// Store dS in euler array
		x_euler_array[atom]=xyz[0];
		y_euler_array[atom]=xyz[1];
		z_euler_array[atom]=xyz[2];

		// Calculate Euler Step
		S_new[0]=S[0]+xyz[0]*mp::dt;
		S_new[1]=S[1]+xyz[1]*mp::dt;
		S_new[2]=S[2]+xyz[2]*mp::dt;

		// Normalise Spin Length
		mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

		S_new[0]=S_new[0]*mod_S;
		S_new[1]=S_new[1]*mod_S;
		S_new[2]=S_new[2]*mod_S;

		//Writing of Spin Values to Storage Array
		x_spin_storage_array[atom]=S_new[0];
		y_spin_storage_array[atom]=S_new[1];
		z_spin_storage_array[atom]=S_new[2];
 	}

	// Copy new spins to spin array
	for(int atom=0;atom<num_atoms;atom++){
		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
		atoms::z_spin_array[atom]=z_spin_storage_array[atom];
	}

	// Recalculate spin dependent fields
	sim::calculate_spin_fields(0,num_atoms);

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

		// Calculate Delta S
		xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
		xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
		xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

		// Store dS in heun array
		x_heun_array[atom]=xyz[0];
		y_heun_array[atom]=xyz[1];
		z_heun_array[atom]=xyz[2];
	}

	// Calculate Heun Step
	for(int atom=0;atom<num_atoms;atom++){
		S_new[0]=x_initial_spin_array[atom]+mp::half_dt*(x_euler_array[atom]+x_heun_array[atom]);
		S_new[1]=y_initial_spin_array[atom]+mp::half_dt*(y_euler_array[atom]+y_heun_array[atom]);
		S_new[2]=z_initial_spin_array[atom]+mp::half_dt*(z_euler_array[atom]+z_heun_array[atom]);

		// Normalise Spin Length
		mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

		S_new[0]=S_new[0]*mod_S;
		S_new[1]=S_new[1]*mod_S;
		S_new[2]=S_new[2]*mod_S;

		// Copy new spins to spin array
		atoms::x_spin_array[atom]=S_new[0];
		atoms::y_spin_array[atom]=S_new[1];
		atoms::z_spin_array[atom]=S_new[2];
	}

	return;
}

} // end of intental namespace

} // end of sim namespace
