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
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

#include <cmath>

int calculate_spin_fields(const int,const int);
int calculate_external_fields(const int,const int);
int set_LLG();

namespace sim{

int LLG_Heun_mpi(){
	//======================================================
	// Subroutine to perform a single LLG integration step
	//======================================================


		//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "LLG_Heun_mpi has been called" << std::endl;}

	using namespace LLG_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGinit();

	//----------------------------------------
	// Local variables for system generation
	//----------------------------------------
	//const int num_atoms = atoms::num_atoms;
	const int pre_comm_si = 0;
	const int pre_comm_ei = vmpi::num_core_atoms;
	const int post_comm_si = vmpi::num_core_atoms;
	const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

	double xyz[3];		/// Local Delta Spin Components
	double S_new[3];	/// New Local Spin Moment
	double mod_S;		/// magnitude of spin moment

		//----------------------------------------
		// Initiate halo swap
		//----------------------------------------
		vmpi::mpi_init_halo_swap();

		//----------------------------------------
		// Store initial spin positions (all)
		//----------------------------------------

		for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
			x_initial_spin_array[atom] = atoms::x_spin_array[atom];
			y_initial_spin_array[atom] = atoms::y_spin_array[atom];
			z_initial_spin_array[atom] = atoms::z_spin_array[atom];
			}

		//----------------------------------------
		// Calculate fields (core)
		//----------------------------------------

		calculate_spin_fields(pre_comm_si,pre_comm_ei);
		calculate_external_fields(pre_comm_si,pre_comm_ei);

		//----------------------------------------
		// Calculate Euler Step (Core)
		//----------------------------------------

		for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){

			const int imaterial=atoms::type_array[atom];
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

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
			S_new[0]=S[0]+xyz[0]*material_parameters::dt;
			S_new[1]=S[1]+xyz[1]*material_parameters::dt;
			S_new[2]=S[2]+xyz[2]*material_parameters::dt;

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

		//----------------------------------------
		// Complete halo swap
		//----------------------------------------
		vmpi::mpi_complete_halo_swap();

		//----------------------------------------
		// Calculate fields (boundary)
		//----------------------------------------

		calculate_spin_fields(post_comm_si,post_comm_ei);
		calculate_external_fields(post_comm_si,post_comm_ei);

		//----------------------------------------
		// Calculate Euler Step (boundary)
		//----------------------------------------

		for(int atom=post_comm_si;atom<post_comm_ei;atom++){

			const int imaterial=atoms::type_array[atom];
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

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
			S_new[0]=S[0]+xyz[0]*material_parameters::dt;
			S_new[1]=S[1]+xyz[1]*material_parameters::dt;
			S_new[2]=S[2]+xyz[2]*material_parameters::dt;

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

		//----------------------------------------
		// Copy new spins to spin array (all)
		//----------------------------------------
		for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
			atoms::x_spin_array[atom]=x_spin_storage_array[atom];
			atoms::y_spin_array[atom]=y_spin_storage_array[atom];
			atoms::z_spin_array[atom]=z_spin_storage_array[atom];
		}

		//------------------------------------------
		// Initiate second halo swap
		//------------------------------------------
		vmpi::mpi_init_halo_swap();

		//------------------------------------------
		// Recalculate spin dependent fields (core)
		//------------------------------------------

		calculate_spin_fields(pre_comm_si,pre_comm_ei);

		//----------------------------------------
		// Calculate Heun Gradients (core)
		//----------------------------------------

		for(int atom=pre_comm_si;atom<pre_comm_ei;atom++){

			const int imaterial=atoms::type_array[atom];;
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

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

		//------------------------------------------
		// Complete second halo swap
		//------------------------------------------
		vmpi::mpi_complete_halo_swap();

		//------------------------------------------
		// Recalculate spin dependent fields (boundary)
		//------------------------------------------

		calculate_spin_fields(post_comm_si,post_comm_ei);

		//----------------------------------------
		// Calculate Heun Gradients (boundary)
		//----------------------------------------

		for(int atom=post_comm_si;atom<post_comm_ei;atom++){

			const int imaterial=atoms::type_array[atom];;
			const double one_oneplusalpha_sq = material_parameters::material[imaterial].one_oneplusalpha_sq;
			const double alpha_oneplusalpha_sq = material_parameters::material[imaterial].alpha_oneplusalpha_sq;

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

		//----------------------------------------
		// Calculate Heun Step
		//----------------------------------------

		for(int atom=pre_comm_si;atom<post_comm_ei;atom++){
			S_new[0]=x_initial_spin_array[atom]+material_parameters::half_dt*(x_euler_array[atom]+x_heun_array[atom]);
			S_new[1]=y_initial_spin_array[atom]+material_parameters::half_dt*(y_euler_array[atom]+y_heun_array[atom]);
			S_new[2]=z_initial_spin_array[atom]+material_parameters::half_dt*(z_euler_array[atom]+z_heun_array[atom]);

			// Normalise Spin Length
			mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

			S_new[0]=S_new[0]*mod_S;
			S_new[1]=S_new[1]*mod_S;
			S_new[2]=S_new[2]*mod_S;

			//----------------------------------------
			// Copy new spins to spin array
			//----------------------------------------
			atoms::x_spin_array[atom]=S_new[0];
			atoms::y_spin_array[atom]=S_new[1];
			atoms::z_spin_array[atom]=S_new[2];
		}

	// Swap timers compute -> wait
	vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

	// Wait for other processors
	vmpi::barrier();

	// Swap timers wait -> compute
	vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

	return EXIT_SUCCESS;
}

} // end of namespace sim
#endif
