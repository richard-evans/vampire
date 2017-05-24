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
#ifdef MPICF
// Vampire Header Files
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <cmath>

int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);
int set_LLG();

namespace sim{

int LLG_Midpoint_mpi(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "LLG_Midpoint_mpi has been called" << std::endl;}

	using namespace LLG_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGinit();

	// Local variables for core / boundary integration
	const int pre_comm_si = 0;
	const int pre_comm_ei = vmpi::num_core_atoms;
	const int post_comm_si = vmpi::num_core_atoms;
	const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

	// Initiate halo swap
	vmpi::mpi_init_halo_swap();

	// Store initial spin positions (all)
	for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
		x_initial_spin_array[atom] = atoms::x_spin_array[atom];
		y_initial_spin_array[atom] = atoms::y_spin_array[atom];
		z_initial_spin_array[atom] = atoms::z_spin_array[atom];
		}

	// Calculate fields (core)
	calculate_spin_fields(pre_comm_si,pre_comm_ei);
	calculate_external_fields(pre_comm_si,pre_comm_ei);

	// Calculate Predictor Step (Core)
	for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){

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


	// Complete halo swap
	vmpi::mpi_complete_halo_swap();

	// Calculate fields (boundary)
	calculate_spin_fields(post_comm_si,post_comm_ei);
	calculate_external_fields(post_comm_si,post_comm_ei);

	// Calculate Predictor Step (boundary)
	for(int atom=post_comm_si;atom<post_comm_ei;atom++){

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

	// Copy new spins to spin array (all)
	for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
		atoms::z_spin_array[atom]=z_spin_storage_array[atom];
	}

	// Initiate second halo swap
	vmpi::mpi_init_halo_swap();

	// Recalculate spin dependent fields (core)
	calculate_spin_fields(pre_comm_si,pre_comm_ei);

	// Calculate Corrector Step (core)
	for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){

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
		x_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[0]*one_minus_beta2FdotF + 2.0*(beta*(F[1]*S[2]-F[2]*S[1]) + F[0]*beta2FdotS));
		y_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[1]*one_minus_beta2FdotF + 2.0*(beta*(F[2]*S[0]-F[0]*S[2]) + F[1]*beta2FdotS));
		z_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[2]*one_minus_beta2FdotF + 2.0*(beta*(F[0]*S[1]-F[1]*S[0]) + F[2]*beta2FdotS));
	}

	// Complete second halo swap
	vmpi::mpi_complete_halo_swap();

	// Recalculate spin dependent fields (boundary)
	calculate_spin_fields(post_comm_si,post_comm_ei);

	// Calculate Corrector Step (boundary)
	for(int atom=post_comm_si;atom<post_comm_ei;atom++){

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
		x_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[0]*one_minus_beta2FdotF + 2.0*(beta*(F[1]*S[2]-F[2]*S[1]) + F[0]*beta2FdotS));
		y_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[1]*one_minus_beta2FdotF + 2.0*(beta*(F[2]*S[0]-F[0]*S[2]) + F[1]*beta2FdotS));
		z_spin_storage_array[atom] = one_o_one_plus_beta2FdotF*(S[2]*one_minus_beta2FdotF + 2.0*(beta*(F[0]*S[1]-F[1]*S[0]) + F[2]*beta2FdotS));
	}

	// Copy new spins to spin array (all)
	for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
		atoms::x_spin_array[atom]=x_spin_storage_array[atom];
		atoms::y_spin_array[atom]=y_spin_storage_array[atom];
		atoms::z_spin_array[atom]=z_spin_storage_array[atom];
	}

	// Wait for other processors
	vmpi::barrier();

	return EXIT_SUCCESS;
}

} // end of namespace sim
#endif
