//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2021. All rights reserved.
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
#include "constants.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmath.hpp"
#include "vio.hpp"

//Internal header file
#include "internal.hpp"

namespace montecarlo{

/// local cmc namespace
namespace cmc{

	std::vector<cmc_material_t> cmc_mask; // class containing properties for constrained atoms
	std::vector<int> mask;                // list of atom <-> mask relations

//------------------------------------------------------------------------------
// Sets up matrices for performing CMC in an arbitrary space
//------------------------------------------------------------------------------
void mask_polar_rot_matrix(std::vector<cmc_material_t>& cmc_mask){

	// loop over all sets of constrained and unconstrained materials
	for( size_t mask=0; mask < cmc_mask.size(); mask++ ){

		const double phi   = cmc_mask[mask].constraint_phi;
		const double theta = cmc_mask[mask].constraint_theta;

		std::vector<double> reference_vector(3);
		std::vector< std::vector<double> > x_rotation_matrix, y_rotation_matrix, z_rotation_matrix, ref_vec;

		const double pi = M_PI;

		//--------------------------------------------------
		// Initialise varibales
		//--------------------------------------------------

		const double ddx = 0.0; //change in angle in degrees
		const double ddy = phi;
		const double ddz = theta;// assumues x = cos(theta)sin(phi), y = sin(theta)sin(phi)

		const double dx = (ddx/360.0)*2.0*pi; //change in angle in radians
		const double dy = (ddy/360.0)*2.0*pi;
		const double dz = (ddz/360.0)*2.0*pi;

		reference_vector[0] = 0.0;
		reference_vector[1] = 0.0;
		reference_vector[2] = 1.0;

		const double sin_x = sin(dx);
		const double cos_x = cos(dx);
		const double sin_y = sin(dy);
		const double cos_y = cos(dy);
		const double sin_z = sin(dz);
		const double cos_z = cos(dz);

		x_rotation_matrix = vmath::set_matrix(3,3);
		y_rotation_matrix = vmath::set_matrix(3,3);
		z_rotation_matrix = vmath::set_matrix(3,3);
		ref_vec = vmath::set_matrix(1,3,reference_vector);

		x_rotation_matrix[0][0] = 1.0;
		x_rotation_matrix[1][0] = 0.0;
		x_rotation_matrix[2][0] = 0.0;
		x_rotation_matrix[0][1] = 0.0;
		x_rotation_matrix[1][1] = cos_x;
		x_rotation_matrix[2][1] = -sin_x;
		x_rotation_matrix[0][2] = 0.0;
		x_rotation_matrix[1][2] = sin_x;
		x_rotation_matrix[2][2] = cos_x;

		y_rotation_matrix[0][0] = cos_y;
		y_rotation_matrix[1][0] = 0.0;
		y_rotation_matrix[2][0] = sin_y;
		y_rotation_matrix[0][1] = 0.0;
		y_rotation_matrix[1][1] = 1.0;
		y_rotation_matrix[2][1] = 0.0;
		y_rotation_matrix[0][2] = -sin_y;
		y_rotation_matrix[1][2] = 0.0;
		y_rotation_matrix[2][2] = cos_y;

		z_rotation_matrix[0][0] = cos_z;
		z_rotation_matrix[1][0] = -sin_z;
		z_rotation_matrix[2][0] = 0.0;
		z_rotation_matrix[0][1] = sin_z;
		z_rotation_matrix[1][1] = cos_z;
		z_rotation_matrix[2][1] = 0.0;
		z_rotation_matrix[0][2] = 0.0;
		z_rotation_matrix[1][2] = 0.0;
		z_rotation_matrix[2][2] = 1.0;

		std::vector< std::vector<double> > polar_matrix    = vmath::matmul(y_rotation_matrix,z_rotation_matrix);
		std::vector< std::vector<double> > polar_matrix_tp = vmath::transpose(polar_matrix);
		std::vector< std::vector<double> > polar_vector    = vmath::matmul(ref_vec,polar_matrix);

		// copy matrices to performance optimised class variables
		for ( int i = 0 ; i < 3 ; i++ ){
			cmc_mask[mask].ppolar_vector[i] = polar_vector[0][i];
			for ( int j = 0 ; j < 3 ; j++ ){
				cmc_mask[mask].ppolar_matrix[i][j]    = polar_matrix[i][j];
				cmc_mask[mask].ppolar_matrix_tp[i][j] = polar_matrix_tp[i][j];
			}
		}

	} // end of loop over materials

} // end of polar rotation initialisation

} // end of cmc namespace

//------------------------------------------------------------------------------
// Creates the rotation matices for the given angle and computes the
// Boltzmann factor for the given temperature.
//------------------------------------------------------------------------------
void initialise_masked_cmc_mc(const int num_sets,                       // number of sets of constrained and unconstrained atoms
										const std::vector<int>& mask,             // unique ID for N sets of atoms with different constraints
										const std::vector<bool>& constrained,     // flag to indicate if atom set with mask index is constrained
										const std::vector<double>& constraints){  // list of 2N vectors listing constraint angles theta (from z) and phi (from x)

	//-----------------------------------------------------
	// initialise masked cmc properties
	//-----------------------------------------------------

	// copy mask to internal storage
	cmc::mask = mask;

	cmc::cmc_mask.resize( num_sets );

	// temporary variable for brevity
	const double pi180 = M_PI / 180.0;

	for( size_t mask_id = 0; mask_id < cmc::cmc_mask.size(); mask_id++ ){
		cmc::cmc_mask[mask_id].constrained      = constrained[   mask_id     ];
		cmc::cmc_mask[mask_id].constraint_phi   = constraints[ 2*mask_id + 0 ];
		cmc::cmc_mask[mask_id].constraint_theta = constraints[ 2*mask_id + 1 ];
		cmc::cmc_mask[mask_id].sx = sin(cmc::cmc_mask[mask_id].constraint_phi * pi180) * cos(cmc::cmc_mask[mask_id].constraint_theta * pi180);
		cmc::cmc_mask[mask_id].sy = sin(cmc::cmc_mask[mask_id].constraint_phi * pi180) * sin(cmc::cmc_mask[mask_id].constraint_theta * pi180);
		cmc::cmc_mask[mask_id].sz = cos(cmc::cmc_mask[mask_id].constraint_phi * pi180);
	}

	// Create rotational matrices for cmc
	cmc::mask_polar_rot_matrix(cmc::cmc_mask);

	// create lists of atoms in each constrained set
	cmc::atom_list.resize(num_sets);
	// create list of matching materials
	for( int atom = 0; atom < atoms::num_atoms; atom++){
		const int mask_id = cmc::mask[atom];
		cmc::atom_list[mask_id].push_back(atom);
	}

	// Initialise all spins along the constraint direction(s)
	for(int atom = 0; atom < atoms::num_atoms; atom++){

		// get mask ID of atom
		const int mask_id = cmc::mask[atom];

		if(cmc::cmc_mask[mask_id].constrained){
			atoms::x_spin_array[atom] = cmc::cmc_mask[mask_id].sx;
			atoms::y_spin_array[atom] = cmc::cmc_mask[mask_id].sy;
			atoms::z_spin_array[atom] = cmc::cmc_mask[mask_id].sz;
		}

	}

	// disable thermal field calculation
	sim::hamiltonian_simulation_flags[3]=0;

	// set initialised flag to true
	cmc::is_initialised = true;
	cmc::masked_cmc = true;

	return;

}

//------------------------------------------------------------------------------
// Chooses nspins random spin pairs from the spin system and attempts a
// Constrained Monte Carlo move on each pair, accepting for either lower
// energy or with a Boltzmann thermal weighting.
//------------------------------------------------------------------------------
void cmc_mc_step_mask(){

	/*double delta_energy1;
	double delta_energy2;
	double delta_energy21;

	double Eold;
	double Enew;*/

   std::vector<double> spin1_initial(3);
	std::vector<double> spin1_final(3);
	double spin2_initial[3];
	double spin2_final[3];

	double spin1_init_mvd[3];
	double spin1_fin_mvd[3];
	double spin2_init_mvd[3];
	double spin2_fin_mvd[3];

	/*double Mz_old;
	double Mz_new;

	double probability;*/

	const double muB = constants::muB;
	const double kB  = constants::kB;

	//---------------------------------------------------------------------------
   // Material dependent temperature rescaling
	//---------------------------------------------------------------------------
   std::vector<double> rescaled_material_kBTBohr(mp::num_materials);
   std::vector<double> sigma_array(mp::num_materials); // range for tuned gaussian random move
   for(int m = 0; m < mp::num_materials; ++m){
      double alpha = mp::material[m].temperature_rescaling_alpha;
      double Tc    = mp::material[m].temperature_rescaling_Tc;
      double rescaled_temperature = sim::temperature < Tc ? Tc*pow(sim::temperature/Tc,alpha) : sim::temperature;
      rescaled_material_kBTBohr[m] = muB/(rescaled_temperature*kB);
      sigma_array[m] = pow(1.0/rescaled_material_kBTBohr[m],0.2)*0.08;
   }

	//---------------------------------------------------------------------------
	// save initial magnetisations
	//---------------------------------------------------------------------------
	for(size_t mask_id = 0; mask_id < cmc::cmc_mask.size(); mask_id++){
		cmc::cmc_mask[mask_id].M_other[0] = 0.0;
		cmc::cmc_mask[mask_id].M_other[1] = 0.0;
		cmc::cmc_mask[mask_id].M_other[2] = 0.0;
	}
	for(int atom = 0; atom < atoms::num_atoms; atom++){
		const int mask_id = cmc::mask[atom];	// get mask ID of atom
		cmc::cmc_mask[mask_id].M_other[0] += atoms::x_spin_array[atom]; // multiplied by polar_vector below
		cmc::cmc_mask[mask_id].M_other[1] += atoms::y_spin_array[atom];
		cmc::cmc_mask[mask_id].M_other[2] += atoms::z_spin_array[atom];
	}

	double statistics_reject = 0.0; // counter for number of rejected moves (for adaptive algorothm)

	//---------------------------------------------------------------------------
	// make a sequence of Monte Carlo moves
	//---------------------------------------------------------------------------
	for (int mcs = 0; mcs < atoms::num_atoms; mcs++){

		// Randomly select spin number 1
		const int atom1 = int( mtrandom::grnd() * atoms::num_atoms );

		// get mask and materiual ID
		const int mat1  = atoms::type_array[atom1];
		const int mask1 = cmc::mask[atom1];

		// set up intenral angle based on material 1
      internal::delta_angle = sigma_array[mat1];

		// check for constrained or unconstrained
		//----------------
		// normal MC
		//----------------
		if(cmc::cmc_mask[mask1].constrained==false){

			// Save old spin position
			spin1_initial[0] = atoms::x_spin_array[atom1];
			spin1_initial[1] = atoms::y_spin_array[atom1];
			spin1_initial[2] = atoms::z_spin_array[atom1];

         // Make Monte Carlo move
         montecarlo::internal::mc_move(spin1_initial, spin1_final);

			// Calculate current energy
			const double Eold = sim::calculate_spin_energy(atom1);

			// Copy new spin position
			atoms::x_spin_array[atom1] = spin1_final[0];
			atoms::y_spin_array[atom1] = spin1_final[1];
			atoms::z_spin_array[atom1] = spin1_final[2];

			// Calculate new energy
			const double Enew = sim::calculate_spin_energy(atom1);

			// Calculate difference in Joules/mu_B
			const double delta_energy1 = (Enew-Eold)*mp::material[mat1].mu_s_SI*1.07828231e23; //1/9.27400915e-24

			// Check for lower energy state and accept unconditionally
			if(delta_energy1 < 0.0) cmc::mc_success += 1.0;

			// Otherwise evaluate probability for move
			else{
				if(exp(-delta_energy1*rescaled_material_kBTBohr[mat1]) >= mtrandom::grnd()) cmc::mc_success += 1.0;
				// If rejected reset spin coordinates and continue
				else{
					atoms::x_spin_array[atom1] = spin1_initial[0];
					atoms::y_spin_array[atom1] = spin1_initial[1];
					atoms::z_spin_array[atom1] = spin1_initial[2];
               cmc::energy_reject += 1.0;
					statistics_reject += 1.0;
				}
			}
		}
		//-----------------------------
		// constrained MC move
		//-----------------------------
		else{

			// Save initial Spin 1
			spin1_initial[0] = atoms::x_spin_array[atom1];
			spin1_initial[1] = atoms::y_spin_array[atom1];
			spin1_initial[2] = atoms::z_spin_array[atom1];

			//spin1_init_mvd = matmul(polar_matrix, spin1_initial)
			spin1_init_mvd[0]=cmc::cmc_mask[mask1].ppolar_matrix[0][0]*spin1_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[0][1]*spin1_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[0][2]*spin1_initial[2];
			spin1_init_mvd[1]=cmc::cmc_mask[mask1].ppolar_matrix[1][0]*spin1_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[1][1]*spin1_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[1][2]*spin1_initial[2];
			spin1_init_mvd[2]=cmc::cmc_mask[mask1].ppolar_matrix[2][0]*spin1_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[2][1]*spin1_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[2][2]*spin1_initial[2];

	      // Make Monte Carlo move
	      montecarlo::internal::mc_move(spin1_initial, spin1_final);

			//spin1_fin_mvd = matmul(polar_matrix, spin1_final)
			spin1_fin_mvd[0]=cmc::cmc_mask[mask1].ppolar_matrix[0][0]*spin1_final[0]+cmc::cmc_mask[mask1].ppolar_matrix[0][1]*spin1_final[1]+cmc::cmc_mask[mask1].ppolar_matrix[0][2]*spin1_final[2];
			spin1_fin_mvd[1]=cmc::cmc_mask[mask1].ppolar_matrix[1][0]*spin1_final[0]+cmc::cmc_mask[mask1].ppolar_matrix[1][1]*spin1_final[1]+cmc::cmc_mask[mask1].ppolar_matrix[1][2]*spin1_final[2];
			spin1_fin_mvd[2]=cmc::cmc_mask[mask1].ppolar_matrix[2][0]*spin1_final[0]+cmc::cmc_mask[mask1].ppolar_matrix[2][1]*spin1_final[1]+cmc::cmc_mask[mask1].ppolar_matrix[2][2]*spin1_final[2];

			// Calculate current energy
			const double Eold = sim::calculate_spin_energy(atom1);

			// Copy new spin position (provisionally accept move)
			atoms::x_spin_array[atom1] = spin1_final[0];
			atoms::y_spin_array[atom1] = spin1_final[1];
			atoms::z_spin_array[atom1] = spin1_final[2];

			// Calculate new energy
			const double Enew = sim::calculate_spin_energy(atom1);

			// Calculate difference in Joules/mu_B
			const double delta_energy1 = (Enew-Eold)*mp::material[mat1].mu_s_SI*1.07828231e23; //1/9.27400915e-24

			// Compute second move

			// Randomly select spin number 2 (i/=j) of same material type
			const int atom2 = cmc::atom_list[mask1][int(mtrandom::grnd()*cmc::atom_list[mask1].size())];
			const int mat2  = atoms::type_array[atom2];

			// Save initial Spin 2
			spin2_initial[0] = atoms::x_spin_array[atom2];
			spin2_initial[1] = atoms::y_spin_array[atom2];
			spin2_initial[2] = atoms::z_spin_array[atom2];

			//spin2_init_mvd = matmul(polar_matrix, spin2_initial)
			spin2_init_mvd[0]=cmc::cmc_mask[mask1].ppolar_matrix[0][0]*spin2_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[0][1]*spin2_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[0][2]*spin2_initial[2];
			spin2_init_mvd[1]=cmc::cmc_mask[mask1].ppolar_matrix[1][0]*spin2_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[1][1]*spin2_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[1][2]*spin2_initial[2];
			spin2_init_mvd[2]=cmc::cmc_mask[mask1].ppolar_matrix[2][0]*spin2_initial[0]+cmc::cmc_mask[mask1].ppolar_matrix[2][1]*spin2_initial[1]+cmc::cmc_mask[mask1].ppolar_matrix[2][2]*spin2_initial[2];

			// Calculate new spin based on constraint Mx=My=0
			spin2_fin_mvd[0] = spin1_init_mvd[0] + spin2_init_mvd[0] - spin1_fin_mvd[0];
			spin2_fin_mvd[1] = spin1_init_mvd[1] + spin2_init_mvd[1] - spin1_fin_mvd[1];

			if(((spin2_fin_mvd[0]*spin2_fin_mvd[0]+spin2_fin_mvd[1]*spin2_fin_mvd[1])<1.0) && (atom1 != atom2)){

				spin2_fin_mvd[2] = vmath::sign(spin2_init_mvd[2])*sqrt(1.0-spin2_fin_mvd[0]*spin2_fin_mvd[0] - spin2_fin_mvd[1]*spin2_fin_mvd[1]);

				//spin2_final = matmul(polar_matrix_tp, spin2_fin_mvd)
				spin2_final[0]=cmc::cmc_mask[mask1].ppolar_matrix_tp[0][0]*spin2_fin_mvd[0]+cmc::cmc_mask[mask1].ppolar_matrix_tp[0][1]*spin2_fin_mvd[1]+cmc::cmc_mask[mask1].ppolar_matrix_tp[0][2]*spin2_fin_mvd[2];
				spin2_final[1]=cmc::cmc_mask[mask1].ppolar_matrix_tp[1][0]*spin2_fin_mvd[0]+cmc::cmc_mask[mask1].ppolar_matrix_tp[1][1]*spin2_fin_mvd[1]+cmc::cmc_mask[mask1].ppolar_matrix_tp[1][2]*spin2_fin_mvd[2];
				spin2_final[2]=cmc::cmc_mask[mask1].ppolar_matrix_tp[2][0]*spin2_fin_mvd[0]+cmc::cmc_mask[mask1].ppolar_matrix_tp[2][1]*spin2_fin_mvd[1]+cmc::cmc_mask[mask1].ppolar_matrix_tp[2][2]*spin2_fin_mvd[2];

				//Calculate Energy Difference 2
				// Calculate current energy
				const double Eold = sim::calculate_spin_energy(atom2);

	         // Copy new spin position (provisionally accept move)
				atoms::x_spin_array[atom2] = spin2_final[0];
				atoms::y_spin_array[atom2] = spin2_final[1];
				atoms::z_spin_array[atom2] = spin2_final[2];

				// Calculate new energy
				const double Enew = sim::calculate_spin_energy(atom2);

	         // Calculate difference in Joules/mu_B
				const double delta_energy2 = (Enew-Eold)*mp::material[mat2].mu_s_SI*1.07828231e23; //1/9.27400915e-24

				// Calculate Delta E for both spins
				const double delta_energy21 = delta_energy1*rescaled_material_kBTBohr[mat1] + delta_energy2*rescaled_material_kBTBohr[mat2];

				// Compute Mz_other, Mz, Mz'
				const double Mz_old = cmc::cmc_mask[mask1].M_other[0] * cmc::cmc_mask[mask1].ppolar_vector[0] +
											 cmc::cmc_mask[mask1].M_other[1] * cmc::cmc_mask[mask1].ppolar_vector[1] +
											 cmc::cmc_mask[mask1].M_other[2] * cmc::cmc_mask[mask1].ppolar_vector[2];

				const double Mz_new = (cmc::cmc_mask[mask1].M_other[0] + spin1_final[0] + spin2_final[0]- spin1_initial[0] - spin2_initial[0])*cmc::cmc_mask[mask1].ppolar_vector[0] +
											 (cmc::cmc_mask[mask1].M_other[1] + spin1_final[1] + spin2_final[1]- spin1_initial[1] - spin2_initial[1])*cmc::cmc_mask[mask1].ppolar_vector[1] +
											 (cmc::cmc_mask[mask1].M_other[2] + spin1_final[2] + spin2_final[2]- spin1_initial[2] - spin2_initial[2])*cmc::cmc_mask[mask1].ppolar_vector[2];

				// Check for lower energy state and accept unconditionally (this allows it to flip sign, why?)
				//if((delta_energy21<0.0) && (Mz_new >= 1e-15) ) continue;

				// If move is favorable then accept
				const double probability = exp(-delta_energy21)*((Mz_new/Mz_old)*(Mz_new/Mz_old))*std::fabs(spin2_init_mvd[2]/spin2_fin_mvd[2]);
				if( ( probability >= mtrandom::grnd() ) && ( Mz_new >= 1e-15 ) ){
					cmc::cmc_mask[mask1].M_other[0] = cmc::cmc_mask[mask1].M_other[0] + spin1_final[0] + spin2_final[0] - spin1_initial[0] - spin2_initial[0];
					cmc::cmc_mask[mask1].M_other[1] = cmc::cmc_mask[mask1].M_other[1] + spin1_final[1] + spin2_final[1] - spin1_initial[1] - spin2_initial[1];
					cmc::cmc_mask[mask1].M_other[2] = cmc::cmc_mask[mask1].M_other[2] + spin1_final[2] + spin2_final[2] - spin1_initial[2] - spin2_initial[2];
					cmc::mc_success += 1.0;
				}
				//if both p1 and p2 not allowed then
				else{
					// reset spin positions
					atoms::x_spin_array[atom1] = spin1_initial[0];
					atoms::y_spin_array[atom1] = spin1_initial[1];
					atoms::z_spin_array[atom1] = spin1_initial[2];

					atoms::x_spin_array[atom2] = spin2_initial[0];
					atoms::y_spin_array[atom2] = spin2_initial[1];
					atoms::z_spin_array[atom2] = spin2_initial[2];

					cmc::energy_reject += 1.0;
					statistics_reject += 1.0;

				}
			}
			// if s2 not on unit sphere
			else{
				atoms::x_spin_array[atom1] = spin1_initial[0];
				atoms::y_spin_array[atom1] = spin1_initial[1];
				atoms::z_spin_array[atom1] = spin1_initial[2];
				cmc::sphere_reject+=1.0;
				statistics_reject += 1.0;
			}
		} // end of cmc move
		cmc::mc_total += 1.0;

	} // end of mc loop

	// calculate new adaptive step sigma angle
	if(montecarlo::algorithm == montecarlo::adaptive){
		const double statistics_moves = atoms::num_atoms;
		const double last_rejection_rate = statistics_reject / statistics_moves;
		const double factor = 0.5 / last_rejection_rate;
		montecarlo::internal::adaptive_sigma *= factor;
		// check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
		if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5) montecarlo::internal::adaptive_sigma = 60.0;
	}

	return;

}

} // End of namespace montecarlo
