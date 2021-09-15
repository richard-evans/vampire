//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmath.hpp"
#include "vio.hpp"

// Internal header files
#include "internal.hpp"

namespace montecarlo{

/// local cmc namespace
namespace cmc{

	// Statistics collection
	double mc_success=0.0;
	double mc_total=0.0;
	double sphere_reject=0.0;
	double energy_reject=0.0;

	bool is_initialised=false;

	// Rotational matrices
	std::vector<std::vector<double> > polar_vector;
	std::vector<std::vector<double> > polar_matrix_tp;
	std::vector<std::vector<double> > polar_matrix;

///
/// @brief Sets up matrices for performing CMC in an arbitrary space
///
/// @param[in]		phi					Magnetisation azimuthal angle
/// @param[in]		theta					Magnetisation rotational angle
/// @param[out]	polar_matrix		Matrix to rotate v to z
/// @param[out]	polar_matrix_tp	Matrix to rotate z to v
/// @param[out]	polar_vector 		Constraint Vector
/// @return			void
///
void polar_rot_matrix(
	double phi,
	double theta,
	std::vector< std::vector<double> > & polar_matrix,
	std::vector< std::vector<double> > & polar_matrix_tp,
	std::vector< std::vector<double> > & polar_vector)
{
	double dx,dy,dz; //change in angle in radians
	double ddx,ddy,ddz; //change in angle in degrees

	double sin_x,sin_y,sin_z,cos_x,cos_y,cos_z;
	std::vector<double> reference_vector(3);

	double pi=3.14159265358979323846264338327;

	//--------------------------------------------------
	// Initialise varibales
	//--------------------------------------------------

	ddx = 0.0;
	ddy = phi;
	ddz = theta;// assumues x = cos(theta)sin(phi), y = sin(theta)sin(phi)

	dx = (ddx/360.0)*2.0*pi;
	dy = (ddy/360.0)*2.0*pi;
	dz = (ddz/360.0)*2.0*pi;

	reference_vector[0] = 0.0;
	reference_vector[1] = 0.0;
	reference_vector[2] = 1.0;

	sin_x = sin(dx);
	cos_x = cos(dx);
	sin_y = sin(dy);
	cos_y = cos(dy);
	sin_z = sin(dz);
	cos_z = cos(dz);

	std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix,ref_vec;

	x_rotation_matrix=vmath::set_matrix(3,3);
	y_rotation_matrix=vmath::set_matrix(3,3);
	z_rotation_matrix=vmath::set_matrix(3,3);
	ref_vec=vmath::set_matrix(1,3,reference_vector);


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

	polar_matrix = vmath::matmul(y_rotation_matrix,z_rotation_matrix);
	polar_matrix_tp = vmath::transpose(polar_matrix);
	polar_vector = vmath::matmul(ref_vec,polar_matrix);

}

/// Function to rotate all spin around the z-axis
void rotate_spins_around_z_axis(double ddz){

	std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix;

	// determine rotational matrices for phi, theta rotation
	vmath::set_rotational_matrix(0.0, 0.0, ddz, x_rotation_matrix,y_rotation_matrix,z_rotation_matrix);

	// loop over all spins and rotate by theta around z
	for(int atom =0;atom<atoms::num_atoms;atom++){

			// Load spin coordinates
			internal::Sold[0]=atoms::x_spin_array[atom];
			internal::Sold[1]=atoms::y_spin_array[atom];
			internal::Sold[2]=atoms::z_spin_array[atom];

			// Calculate new spin positions
			internal::Snew = vmath::matmul(internal::Sold,z_rotation_matrix);

			// Set new spin positions
			atoms::x_spin_array[atom]=internal::Snew[0];
			atoms::y_spin_array[atom]=internal::Snew[1];
			atoms::z_spin_array[atom]=internal::Snew[2];
		}

	return;
}

/// Function to rotate all spin around the x-axis
void rotate_spins_around_x_axis(double ddx){

	std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix;

	// determine rotational matrices for phi, theta rotation
	vmath::set_rotational_matrix(ddx, 0.0, 0.0, x_rotation_matrix,y_rotation_matrix,z_rotation_matrix);

	//vmath::print_matrix(x_rotation_matrix);

	// loop over all spins and rotate by phi around x
	for(int atom =0;atom<atoms::num_atoms;atom++){

		// Load spin coordinates
		internal::Sold[0]=atoms::x_spin_array[atom];
		internal::Sold[1]=atoms::y_spin_array[atom];
		internal::Sold[2]=atoms::z_spin_array[atom];

		// Calculate new spin positions
		internal::Snew = vmath::matmul(internal::Sold,x_rotation_matrix);

		// Set new spin positions
		atoms::x_spin_array[atom]=internal::Snew[0];
		atoms::y_spin_array[atom]=internal::Snew[1];
		atoms::z_spin_array[atom]=internal::Snew[2];
	}

	return;
}

} // end of cmc namespace

//------------------------------------------------------------------------------
// Creates the rotation matices for the given angle and computes the
// Boltzmann factor for the given temperature.
//------------------------------------------------------------------------------
void CMCinit(){

	// Check for calling of function
	if(err::check==true) std::cout << "montecarlo::CMCinit has been called" << std::endl;

	// Create rotational matrices for cmc
	cmc::polar_rot_matrix(sim::constraint_phi,sim::constraint_theta, cmc::polar_matrix, cmc::polar_matrix_tp, cmc::polar_vector);

	// Check for rotational update
	if((sim::constraint_rotation==false || (sim::constraint_theta_changed==false && sim::constraint_phi_changed==false)) && sim::checkpoint_loaded_flag==false){

		// Output message showing constraint direction re-initialisation
		zlog << zTs() << "Initialising spins to new constraint direction (phi, theta) " <<  sim::constraint_phi << " , " << sim::constraint_theta << std::endl;

		// Initialise all spins along the constraint direction.
		double sx=sin(sim::constraint_phi*M_PI/180.0)*cos(sim::constraint_theta*M_PI/180.0);
		double sy=sin(sim::constraint_phi*M_PI/180.0)*sin(sim::constraint_theta*M_PI/180.0);
		double sz=cos(sim::constraint_phi*M_PI/180.0);

		for(int atom =0;atom<atoms::num_atoms;atom++){
			atoms::x_spin_array[atom]=sx;
			atoms::y_spin_array[atom]=sy;
			atoms::z_spin_array[atom]=sz;
		}
	}
   // If the simulation is restarted from within the temperature loop, then the spin configurations must not be reinitialised
   else if(sim::checkpoint_loaded_flag==true){
      zlog << zTs() << "Not re-initialising spins to constraint directions since still in temperature loop" << std::endl;
   }
	else{

		// Output message showing constraint direction re-initialisation
		zlog << zTs() << "Initialising spins by rotation from last constraint direction to new constraint direction (phi, theta) "
		<<  sim::constraint_phi << " , " << sim::constraint_theta << std::endl;

		// Rotate spins from old to new constraint direction

		// Determine angles to rotate spins by
		double theta_old = sim::constraint_theta;
		double theta_new = sim::constraint_theta;

		double phi_old = sim::constraint_phi;
		double phi_new = sim::constraint_phi;

		// note - sim::constraint_phi, sim::constraint_theta are already at new angle
		if(sim::constraint_theta_changed) theta_old = sim::constraint_theta - sim::constraint_theta_delta;
		if(sim::constraint_phi_changed) phi_old     = sim::constraint_phi   - sim::constraint_phi_delta;

		// Rotate all spins from theta_old to theta = 0 (reference direction along x)
		cmc::rotate_spins_around_z_axis(-theta_old);

		// Rotate all spins from phi_old to phi_new
		cmc::rotate_spins_around_x_axis(phi_new-phi_old);

		// Rotate all spins from theta = 0 to theta = theta_new
		cmc::rotate_spins_around_z_axis(theta_new);

		// reset rotation flags
		sim::constraint_theta_changed = false;
		sim::constraint_phi_changed   = false;

	}

	// disable thermal field calculation
	sim::hamiltonian_simulation_flags[3]=0;

	// set initialised flag to true
	cmc::is_initialised=true;

	return;
}

//------------------------------------------------------------------------------
// Chooses nspins random spin pairs from the spin system and attempts a
// Constrained Monte Carlo move on each pair, accepting for either lower
// energy or with a Boltzmann thermal weighting.
//------------------------------------------------------------------------------
int cmc_step(){

	// Check for calling of function
	if(err::check==true) std::cout << "montecarlo::ConstrainedMonteCarlo has been called" << std::endl;

	// check for cmc initialisation
	if(cmc::is_initialised==false) CMCinit();

	int atom_number1;
	int atom_number2;
	int imat1;
	int imat2;

	double delta_energy1;
	double delta_energy2;
	double delta_energy21;

	double Eold;
	double Enew;

	std::vector<double> spin1_initial(3);
	std::vector<double> spin1_final(3);
	double spin2_initial[3];
	double spin2_final[3];

	double spin1_init_mvd[3];
	double spin1_fin_mvd[3];
	double spin2_init_mvd[3];
	double spin2_fin_mvd[3];

	double M_other[3];
	double Mz_old;
	double Mz_new;

	double probability;

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

	// copy matrices for speed
	double ppolar_vector[3];
	double ppolar_matrix[3][3];
	double ppolar_matrix_tp[3][3];

	for (int i=0;i<3;i++){
		ppolar_vector[i]=cmc::polar_vector[0][i];
		for (int j=0;j<3;j++){
			ppolar_matrix[i][j]=cmc::polar_matrix[i][j];
			ppolar_matrix_tp[i][j]=cmc::polar_matrix_tp[i][j];
		}
	}

	//sigma = ((kbT_bohr)**0.2)*0.08_dp
	//inv_kbT_bohr = 1.0_dp/kbT_bohr

	M_other[0] = 0.0;
	M_other[1] = 0.0;
	M_other[2] = 0.0;

	for(int atom=0;atom<atoms::num_atoms;atom++){
		const double mu = atoms::m_spin_array[atom];
		M_other[0] = M_other[0] + atoms::x_spin_array[atom] * mu; //multiplied by polar_vector below
		M_other[1] = M_other[1] + atoms::y_spin_array[atom] * mu;
		M_other[2] = M_other[2] + atoms::z_spin_array[atom] * mu;
	}

	double statistics_reject = 0.0; // counter for number of rejected moves (for adaptive algorothm)

	for (int mcs=0;mcs<atoms::num_atoms;mcs++){
		// Randomly select spin number 1
		atom_number1 = int(mtrandom::grnd()*atoms::num_atoms);
		imat1=atoms::type_array[atom_number1];

		// get spin moment for atom 1
		const double mu1 = atoms::m_spin_array[atom_number1];

      internal::delta_angle=sigma_array[imat1];

		// Save initial Spin 1
		spin1_initial[0] = atoms::x_spin_array[atom_number1];
		spin1_initial[1] = atoms::y_spin_array[atom_number1];
		spin1_initial[2] = atoms::z_spin_array[atom_number1];

		//spin1_init_mvd = matmul(polar_matrix, spin1_initial)
		spin1_init_mvd[0]=ppolar_matrix[0][0]*spin1_initial[0]+ppolar_matrix[0][1]*spin1_initial[1]+ppolar_matrix[0][2]*spin1_initial[2];
		spin1_init_mvd[1]=ppolar_matrix[1][0]*spin1_initial[0]+ppolar_matrix[1][1]*spin1_initial[1]+ppolar_matrix[1][2]*spin1_initial[2];
		spin1_init_mvd[2]=ppolar_matrix[2][0]*spin1_initial[0]+ppolar_matrix[2][1]*spin1_initial[1]+ppolar_matrix[2][2]*spin1_initial[2];

      // Make Monte Carlo move
      montecarlo::internal::mc_move(spin1_initial,spin1_final);

		//spin1_fin_mvd = matmul(polar_matrix, spin1_final)
		spin1_fin_mvd[0]=ppolar_matrix[0][0]*spin1_final[0]+ppolar_matrix[0][1]*spin1_final[1]+ppolar_matrix[0][2]*spin1_final[2];
		spin1_fin_mvd[1]=ppolar_matrix[1][0]*spin1_final[0]+ppolar_matrix[1][1]*spin1_final[1]+ppolar_matrix[1][2]*spin1_final[2];
		spin1_fin_mvd[2]=ppolar_matrix[2][0]*spin1_final[0]+ppolar_matrix[2][1]*spin1_final[1]+ppolar_matrix[2][2]*spin1_final[2];

		// Calculate Energy Difference 1
		//call calc_one_spin_energy(delta_energy1,spin1_final,atom_number1)

		// Calculate current energy
		Eold = sim::calculate_spin_energy(atom_number1);

		// Copy new spin position (provisionally accept move)
		atoms::x_spin_array[atom_number1] = spin1_final[0];
		atoms::y_spin_array[atom_number1] = spin1_final[1];
		atoms::z_spin_array[atom_number1] = spin1_final[2];

		// Calculate new energy
		Enew = sim::calculate_spin_energy(atom_number1);

		// Calculate difference in Joules/mu_B
		delta_energy1 = (Enew-Eold)*mp::material[imat1].mu_s_SI*1.07828231e23; //1/9.27400915e-24

		// Compute second move

		// Randomly select spin number 2 (i/=j)
		atom_number2 = int(mtrandom::grnd()*atoms::num_atoms);
		imat2=atoms::type_array[atom_number2];
      // get spin moment for atom 2
      const double mu2 = atoms::m_spin_array[atom_number2];

		// Save initial Spin 2
		spin2_initial[0] = atoms::x_spin_array[atom_number2];
		spin2_initial[1] = atoms::y_spin_array[atom_number2];
		spin2_initial[2] = atoms::z_spin_array[atom_number2];

		//spin2_init_mvd = matmul(polar_matrix, spin2_initial)
		spin2_init_mvd[0]=ppolar_matrix[0][0]*spin2_initial[0]+ppolar_matrix[0][1]*spin2_initial[1]+ppolar_matrix[0][2]*spin2_initial[2];
		spin2_init_mvd[1]=ppolar_matrix[1][0]*spin2_initial[0]+ppolar_matrix[1][1]*spin2_initial[1]+ppolar_matrix[1][2]*spin2_initial[2];
		spin2_init_mvd[2]=ppolar_matrix[2][0]*spin2_initial[0]+ppolar_matrix[2][1]*spin2_initial[1]+ppolar_matrix[2][2]*spin2_initial[2];

		// Calculate new spin based on constraint Mx=My=0
		const double mom_ratio = (mu1/mu2);
		spin2_fin_mvd[0] = mom_ratio*(spin1_init_mvd[0]-spin1_fin_mvd[0]) + spin2_init_mvd[0];
		spin2_fin_mvd[1] = mom_ratio*(spin1_init_mvd[1]-spin1_fin_mvd[1]) + spin2_init_mvd[1];

		if(((spin2_fin_mvd[0]*spin2_fin_mvd[0]+spin2_fin_mvd[1]*spin2_fin_mvd[1])<1.0) && (atom_number1 != atom_number2)){

			spin2_fin_mvd[2] = vmath::sign(spin2_init_mvd[2])*sqrt(1.0-spin2_fin_mvd[0]*spin2_fin_mvd[0] - spin2_fin_mvd[1]*spin2_fin_mvd[1]);

			//spin2_final = matmul(polar_matrix_tp, spin2_fin_mvd)
			spin2_final[0]=ppolar_matrix_tp[0][0]*spin2_fin_mvd[0]+ppolar_matrix_tp[0][1]*spin2_fin_mvd[1]+ppolar_matrix_tp[0][2]*spin2_fin_mvd[2];
			spin2_final[1]=ppolar_matrix_tp[1][0]*spin2_fin_mvd[0]+ppolar_matrix_tp[1][1]*spin2_fin_mvd[1]+ppolar_matrix_tp[1][2]*spin2_fin_mvd[2];
			spin2_final[2]=ppolar_matrix_tp[2][0]*spin2_fin_mvd[0]+ppolar_matrix_tp[2][1]*spin2_fin_mvd[1]+ppolar_matrix_tp[2][2]*spin2_fin_mvd[2];

			// Automatically accept move for Spin1 - now before if ^^
			//atomic_spin_array(:,atom_number1) = spin1_final(:)

			//Calculate Energy Difference 2
			// Calculate current energy
			Eold = sim::calculate_spin_energy(atom_number2);

			// Copy new spin position (provisionally accept move)
			atoms::x_spin_array[atom_number2] = spin2_final[0];
			atoms::y_spin_array[atom_number2] = spin2_final[1];
			atoms::z_spin_array[atom_number2] = spin2_final[2];

			// Calculate new energy
			Enew = sim::calculate_spin_energy(atom_number2);

			// Calculate difference in Joules/mu_B
			delta_energy2 = (Enew-Eold)*mp::material[imat2].mu_s_SI*1.07828231e23; //1/9.27400915e-24

			// Calculate Delta E for both spins
			delta_energy21 = delta_energy1*rescaled_material_kBTBohr[imat1] + delta_energy2*rescaled_material_kBTBohr[imat2];

			// Compute Mz_other, Mz, Mz'
			Mz_old = M_other[0]*ppolar_vector[0] + M_other[1]*ppolar_vector[1] + M_other[2]*ppolar_vector[2];

			Mz_new = (M_other[0] + mu1*spin1_final[0] + mu2*spin2_final[0]- mu1*spin1_initial[0] - mu2*spin2_initial[0])*ppolar_vector[0] +
						(M_other[1] + mu1*spin1_final[1] + mu2*spin2_final[1]- mu1*spin1_initial[1] - mu2*spin2_initial[1])*ppolar_vector[1] +
						(M_other[2] + mu1*spin1_final[2] + mu2*spin2_final[2]- mu1*spin1_initial[2] - mu2*spin2_initial[2])*ppolar_vector[2];

			// Check for lower energy state and accept unconditionally
			//if((delta_energy21<0.0) && (Mz_new>0.0)) continue;

			// Otherwise evaluate probability for move
			//else{
				// If move is favorable then accept
				probability = exp(-delta_energy21)*((Mz_new/Mz_old)*(Mz_new/Mz_old))*std::fabs(spin2_init_mvd[2]/spin2_fin_mvd[2]);
				if((probability>=mtrandom::grnd()) && (Mz_new>0.0) ){
					M_other[0] = M_other[0] + mu1*spin1_final[0] + mu2*spin2_final[0] - mu1*spin1_initial[0] - mu2*spin2_initial[0];
					M_other[1] = M_other[1] + mu1*spin1_final[1] + mu2*spin2_final[1] - mu1*spin1_initial[1] - mu2*spin2_initial[1];
					M_other[2] = M_other[2] + mu1*spin1_final[2] + mu2*spin2_final[2] - mu1*spin1_initial[2] - mu2*spin2_initial[2];
					cmc::mc_success += 1.0;
				}
				//if both p1 and p2 not allowed then
				else{
					// reset spin positions
					atoms::x_spin_array[atom_number1] = spin1_initial[0];
					atoms::y_spin_array[atom_number1] = spin1_initial[1];
					atoms::z_spin_array[atom_number1] = spin1_initial[2];

					atoms::x_spin_array[atom_number2] = spin2_initial[0];
					atoms::y_spin_array[atom_number2] = spin2_initial[1];
					atoms::z_spin_array[atom_number2] = spin2_initial[2];

					cmc::energy_reject += 1.0;
					statistics_reject += 1.0;
				}
			//}
		}
		// if s2 not on unit sphere
		else{
			atoms::x_spin_array[atom_number1] = spin1_initial[0];
			atoms::y_spin_array[atom_number1] = spin1_initial[1];
			atoms::z_spin_array[atom_number1] = spin1_initial[2];
			cmc::sphere_reject+=1.0;
			statistics_reject += 1.0;
		}

		cmc::mc_total += 1.0;
	}

	// calculate new adaptive step sigma angle
	if(montecarlo::algorithm == montecarlo::adaptive){
		const double statistics_moves = atoms::num_atoms;
		const double last_rejection_rate = statistics_reject / statistics_moves;
		const double factor = 0.5 / last_rejection_rate;
		montecarlo::internal::adaptive_sigma *= factor;
		// check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
		if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5) montecarlo::internal::adaptive_sigma = 60.0;
	}

	return EXIT_SUCCESS;
}

} // End of namespace montecarlo
