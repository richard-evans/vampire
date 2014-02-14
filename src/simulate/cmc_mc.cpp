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
///	@file  cmc.cpp
///  	@brief  Constrained Monte Carlo / Monte Carlo hybrid integrator
///
/// Additional to CMC, this version combines constrained and free integrators
/// for different materials
///
///	@author Richard Evans, richard.evans@york.ac.uk
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
///  @internal
///    Created  02/11/2011
///   Revision  1.0
///  Copyright  Copyright (c) 2011, Richard Evans
///
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
#include "vmath.hpp"
#include "vio.hpp"

/// local cmc namespace
namespace cmc{
	
	std::vector<std::vector< int > > atom_list;
	std::vector<cmc_material_t> cmc_mat;
	int active_material=0;

///
/// @brief Sets up matrices for performing CMC in an arbitrary space
///
void mat_polar_rot_matrix()
{
	
	for(int mat=0;mat<mp::num_materials;mat++){
		
	const double phi=cmc::cmc_mat[mat].constraint_phi;
	const double theta=cmc::cmc_mat[mat].constraint_theta;
	
	std::vector< std::vector<double> > polar_matrix; 
	std::vector< std::vector<double> > polar_matrix_tp; 
	std::vector< std::vector<double> > polar_vector;
	
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

	// copy matrices to performance optimised class variables
	for (int i=0;i<3;i++){
		cmc::cmc_mat[mat].ppolar_vector[i]=polar_vector[0][i];
		for (int j=0;j<3;j++){
			cmc::cmc_mat[mat].ppolar_matrix[i][j]=polar_matrix[i][j];
			cmc::cmc_mat[mat].ppolar_matrix_tp[i][j]=polar_matrix_tp[i][j];
		}
	}
	
	} // end of loop over materials
	
} // end of polar rotation initialisation

/// Function to rotate all spin around the z-axis
void rotate_material_spins_around_z_axis(double ddz, int material){

	std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix;

	// determine rotational matrices for phi, theta rotation
	vmath::set_rotational_matrix(0.0, 0.0, ddz, x_rotation_matrix,y_rotation_matrix,z_rotation_matrix);

	// loop over all spins and rotate by theta around z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		int mat=atoms::type_array[atom];
		if(mat==material){
			std::vector<double> Sold(3), Snew(3); // Vectors to hold spins

			// Load spin coordinates
			Sold[0]=atoms::x_spin_array[atom];
			Sold[1]=atoms::y_spin_array[atom];
			Sold[2]=atoms::z_spin_array[atom];

			// Calculate new spin positions
			Snew = vmath::matmul(Sold,z_rotation_matrix);

			// Set new spin positions
			atoms::x_spin_array[atom]=Snew[0];
			atoms::y_spin_array[atom]=Snew[1];
			atoms::z_spin_array[atom]=Snew[2];
		}
	}

	return;
}

/// Function to rotate all spin around the x-axis
void rotate_material_spins_around_x_axis(double ddx, int material){

	std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix;

	// determine rotational matrices for phi, theta rotation
	vmath::set_rotational_matrix(ddx, 0.0, 0.0, x_rotation_matrix,y_rotation_matrix,z_rotation_matrix);

	// loop over all spins and rotate by phi around x
	for(int atom =0;atom<atoms::num_atoms;atom++){
		int mat=atoms::type_array[atom];
		if(mat==material){
			std::vector<double> Sold(3), Snew(3); // Vectors to hold spins

			// Load spin coordinates
			Sold[0]=atoms::x_spin_array[atom];
			Sold[1]=atoms::y_spin_array[atom];
			Sold[2]=atoms::z_spin_array[atom];

			// Calculate new spin positions
			Snew = vmath::matmul(Sold,x_rotation_matrix);

			// Set new spin positions
			atoms::x_spin_array[atom]=Snew[0];
			atoms::y_spin_array[atom]=Snew[1];
			atoms::z_spin_array[atom]=Snew[2];
		}
	}

	return;
}

} // end of cmc namespace

namespace sim{
///
/// @brief        Initialise Constrained Monte Carlo module
///
/// Creates the rotation matices for the given angle and computes the
/// Boltzmann factor for the given temperature.
///
/// @return       void
///
void CMCMCinit(){

	// Check for calling of function
	if(err::check==true) std::cout << "sim::CMCMCinit has been called" << std::endl;
	
	// Create rotational matrices for cmc
	cmc::mat_polar_rot_matrix();

	cmc::atom_list.resize(mp::num_materials);
	// create list of matching materials
	for(int atom=0;atom<atoms::num_atoms;atom++){
		int mat=atoms::type_array[atom];
		cmc::atom_list[mat].push_back(atom);
	}
	
		// Check for rotational update
	if(sim::constraint_rotation==false || (sim::constraint_theta_changed==false && sim::constraint_phi_changed==false)){

		// Output message showing constraint direction re-initialisation
		zlog << zTs() << "Initialising spins in all materials to new constraint directions." << std::endl;

	// Initialise all spins along the constraint direction(s).
	for(int atom =0;atom<atoms::num_atoms;atom++){
		int imat=atoms::type_array[atom];
		if(mp::material[imat].constrained==true){
			double sx=sin(cmc::cmc_mat[imat].constraint_phi*M_PI/180.0)*cos(cmc::cmc_mat[imat].constraint_theta*M_PI/180.0);
			double sy=sin(cmc::cmc_mat[imat].constraint_phi*M_PI/180.0)*sin(cmc::cmc_mat[imat].constraint_theta*M_PI/180.0);
			double sz=cos(cmc::cmc_mat[imat].constraint_phi*M_PI/180.0);
			atoms::x_spin_array[atom]=sx;
			atoms::y_spin_array[atom]=sy;
			atoms::z_spin_array[atom]=sz;
		}
	}
	
	}
	else{

		// Output message showing constraint direction re-initialisation
		zlog << zTs() << "Initialising spins in material " << cmc::active_material << " by rotation to new constraint direction (phi, theta) " 
		<<  cmc::cmc_mat[cmc::active_material].constraint_phi << " , " << cmc::cmc_mat[cmc::active_material].constraint_theta << std::endl;

		// Rotate spins from old to new constraint direction

		// Determine angles to rotate spins by
		double theta_old = cmc::cmc_mat[cmc::active_material].constraint_theta; 
		double theta_new = cmc::cmc_mat[cmc::active_material].constraint_theta;

		double phi_old = cmc::cmc_mat[cmc::active_material].constraint_phi;
		double phi_new = cmc::cmc_mat[cmc::active_material].constraint_phi;

		// note - sim::constraint_phi, sim::constraint_theta are already at new angle
		if(sim::constraint_theta_changed) theta_old = cmc::cmc_mat[cmc::active_material].constraint_theta - cmc::cmc_mat[cmc::active_material].constraint_theta_delta;
		if(sim::constraint_phi_changed) phi_old     = cmc::cmc_mat[cmc::active_material].constraint_phi - cmc::cmc_mat[cmc::active_material].constraint_phi_delta;
		
		// Rotate all spins in active material from theta_old to theta = 0 (reference direction along x)
		cmc::rotate_material_spins_around_z_axis(-theta_old, cmc::active_material);

		// Rotate all spins in active material from phi_old to phi_new
		cmc::rotate_material_spins_around_x_axis(phi_new-phi_old, cmc::active_material);

		// Rotate all spins in active material from theta = 0 to theta = theta_new
		cmc::rotate_material_spins_around_z_axis(theta_new, cmc::active_material);

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

///
/// @brief      Runs the Constrained Monte Carlo algorithm
///
/// Chooses nspins random spin pairs from the spin system and attempts a
/// Constrained Monte Carlo move on each pair, accepting for either lower
/// energy or with a Boltzmann thermal weighting.
///
/// @return     void
///
int ConstrainedMonteCarloMonteCarlo(){
	
	// Check for calling of function
	if(err::check==true) std::cout << "sim::ConstrainedMonteCarlo has been called" << std::endl;

	// check for cmc initialisation
	if(cmc::is_initialised==false) CMCMCinit();
	
	int atom_number1;
	int atom_number2;
	int imat1;
	int imat2;

	double delta_energy1;
	double delta_energy2;
	double delta_energy21;

	double Eold;
	double Enew;

   std::valarray<double> spin1_initial(3);
	std::valarray<double> spin1_final(3);
	double spin2_initial[3];
	double spin2_final[3];

	double spin1_init_mvd[3];
	double spin1_fin_mvd[3];
	double spin2_init_mvd[3];
	double spin2_fin_mvd[3];

	double M_other[3];
	double Mz_old;
	double Mz_new;

	double sqrt_ran;
	double probability;
	
	const double kBTBohr = 9.27400915e-24/(sim::temperature*1.3806503e-23);

   // Calculate range for move
   sim::mc_delta_angle=pow(1.0/kBTBohr,0.2)*0.08;
	
	const int AtomExchangeType=atoms::exchange_type; // Cast as constant and pass to energy calculation for speed
	
	// save initial magnetisations
	for(int mat=0;mat<mp::num_materials;mat++){
	cmc::cmc_mat[mat].M_other[0] = 0.0;
	cmc::cmc_mat[mat].M_other[1] = 0.0;
	cmc::cmc_mat[mat].M_other[2] = 0.0;
	}
	for(int atom=0;atom<atoms::num_atoms;atom++){
		int mat=atoms::type_array[atom];
		cmc::cmc_mat[mat].M_other[0] += atoms::x_spin_array[atom]; //multiplied by polar_vector below
		cmc::cmc_mat[mat].M_other[1] += atoms::y_spin_array[atom];
		cmc::cmc_mat[mat].M_other[2] += atoms::z_spin_array[atom];
	}

	// make a sequence of Monte Carlo moves
	for (int mcs=0;mcs<atoms::num_atoms;mcs++){ 
		// Randomly select spin number 1
		atom_number1 = int(mtrandom::grnd()*atoms::num_atoms);
		imat1=atoms::type_array[atom_number1];
		
		// check for constrained or unconstrained
		if(mp::material[imat1].constrained==false){
			// normal MC

			// Save old spin position
			spin1_initial[0] = atoms::x_spin_array[atom_number1];
			spin1_initial[1] = atoms::y_spin_array[atom_number1];
			spin1_initial[2] = atoms::z_spin_array[atom_number1];

         // Make Monte Carlo move
         spin1_final=sim::mc_move(spin1_initial);

			// Calculate current energy
			Eold = sim::calculate_spin_energy(atom_number1, AtomExchangeType);
			
			// Copy new spin position
			atoms::x_spin_array[atom_number1] = spin1_final[0];
			atoms::y_spin_array[atom_number1] = spin1_final[1];
			atoms::z_spin_array[atom_number1] = spin1_final[2];

			// Calculate new energy
			Enew = sim::calculate_spin_energy(atom_number1, AtomExchangeType);
			
			// Calculate difference in Joules/mu_B
			delta_energy1 = (Enew-Eold)*mp::material[imat1].mu_s_SI*1.07828231e23; //1/9.27400915e-24
			
			// Check for lower energy state and accept unconditionally
			if(delta_energy1<0){
            cmc::mc_success += 1.0;
         }
			// Otherwise evaluate probability for move
			else{
				if(exp(-delta_energy1*kBTBohr) >= mtrandom::grnd()){
               cmc::mc_success += 1.0;
            }
				// If rejected reset spin coordinates and continue
				else{
					atoms::x_spin_array[atom_number1] = spin1_initial[0];
					atoms::y_spin_array[atom_number1] = spin1_initial[1];
					atoms::z_spin_array[atom_number1] = spin1_initial[2];
               cmc::energy_reject += 1.0;
				}
			}
		}
		else{
		// constrained MC move
		const int imat=imat1;
		
		// Save initial Spin 1
		spin1_initial[0] = atoms::x_spin_array[atom_number1];
		spin1_initial[1] = atoms::y_spin_array[atom_number1];
		spin1_initial[2] = atoms::z_spin_array[atom_number1];
		
		//spin1_init_mvd = matmul(polar_matrix, spin1_initial)
		spin1_init_mvd[0]=cmc::cmc_mat[imat].ppolar_matrix[0][0]*spin1_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[0][1]*spin1_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[0][2]*spin1_initial[2];
		spin1_init_mvd[1]=cmc::cmc_mat[imat].ppolar_matrix[1][0]*spin1_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[1][1]*spin1_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[1][2]*spin1_initial[2];
		spin1_init_mvd[2]=cmc::cmc_mat[imat].ppolar_matrix[2][0]*spin1_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[2][1]*spin1_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[2][2]*spin1_initial[2];

      // Make Monte Carlo move
      spin1_final=sim::mc_move(spin1_initial);

		//spin1_fin_mvd = matmul(polar_matrix, spin1_final)
		spin1_fin_mvd[0]=cmc::cmc_mat[imat].ppolar_matrix[0][0]*spin1_final[0]+cmc::cmc_mat[imat].ppolar_matrix[0][1]*spin1_final[1]+cmc::cmc_mat[imat].ppolar_matrix[0][2]*spin1_final[2];
		spin1_fin_mvd[1]=cmc::cmc_mat[imat].ppolar_matrix[1][0]*spin1_final[0]+cmc::cmc_mat[imat].ppolar_matrix[1][1]*spin1_final[1]+cmc::cmc_mat[imat].ppolar_matrix[1][2]*spin1_final[2];
		spin1_fin_mvd[2]=cmc::cmc_mat[imat].ppolar_matrix[2][0]*spin1_final[0]+cmc::cmc_mat[imat].ppolar_matrix[2][1]*spin1_final[1]+cmc::cmc_mat[imat].ppolar_matrix[2][2]*spin1_final[2];

		// Calculate current energy
		Eold = sim::calculate_spin_energy(atom_number1, AtomExchangeType);
			
		// Copy new spin position (provisionally accept move)
		atoms::x_spin_array[atom_number1] = spin1_final[0];
		atoms::y_spin_array[atom_number1] = spin1_final[1];
		atoms::z_spin_array[atom_number1] = spin1_final[2];

		// Calculate new energy
		Enew = sim::calculate_spin_energy(atom_number1, AtomExchangeType);
			
		// Calculate difference in Joules/mu_B
		delta_energy1 = (Enew-Eold)*mp::material[imat1].mu_s_SI*1.07828231e23; //1/9.27400915e-24

		// Compute second move

		// Randomly select spin number 2 (i/=j) of same material type
		atom_number2 = cmc::atom_list[imat1][int(mtrandom::grnd()*cmc::atom_list[imat1].size())];
		imat2=atoms::type_array[atom_number2];
		if(imat1!=imat2){
			terminaltextcolor(RED);
			std::cerr << "Error in MC/CMC integration! - atoms pairs are not from same material!" << std::endl;
			terminaltextcolor(WHITE);
			err::vexit();
		}
		
		// Save initial Spin 2
		spin2_initial[0] = atoms::x_spin_array[atom_number2];
		spin2_initial[1] = atoms::y_spin_array[atom_number2];
		spin2_initial[2] = atoms::z_spin_array[atom_number2];
		
		//spin2_init_mvd = matmul(polar_matrix, spin2_initial)
		spin2_init_mvd[0]=cmc::cmc_mat[imat].ppolar_matrix[0][0]*spin2_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[0][1]*spin2_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[0][2]*spin2_initial[2];
		spin2_init_mvd[1]=cmc::cmc_mat[imat].ppolar_matrix[1][0]*spin2_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[1][1]*spin2_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[1][2]*spin2_initial[2];
		spin2_init_mvd[2]=cmc::cmc_mat[imat].ppolar_matrix[2][0]*spin2_initial[0]+cmc::cmc_mat[imat].ppolar_matrix[2][1]*spin2_initial[1]+cmc::cmc_mat[imat].ppolar_matrix[2][2]*spin2_initial[2];

		// Calculate new spin based on constraint Mx=My=0
		spin2_fin_mvd[0] = spin1_init_mvd[0]+spin2_init_mvd[0]-spin1_fin_mvd[0];
		spin2_fin_mvd[1] = spin1_init_mvd[1]+spin2_init_mvd[1]-spin1_fin_mvd[1];

		if(((spin2_fin_mvd[0]*spin2_fin_mvd[0]+spin2_fin_mvd[1]*spin2_fin_mvd[1])<1.0) && (atom_number1 != atom_number2)){ 

			spin2_fin_mvd[2] = vmath::sign(spin2_init_mvd[2])*sqrt(1.0-spin2_fin_mvd[0]*spin2_fin_mvd[0] - spin2_fin_mvd[1]*spin2_fin_mvd[1]);

			//spin2_final = matmul(polar_matrix_tp, spin2_fin_mvd)
			spin2_final[0]=cmc::cmc_mat[imat].ppolar_matrix_tp[0][0]*spin2_fin_mvd[0]+cmc::cmc_mat[imat].ppolar_matrix_tp[0][1]*spin2_fin_mvd[1]+cmc::cmc_mat[imat].ppolar_matrix_tp[0][2]*spin2_fin_mvd[2];
			spin2_final[1]=cmc::cmc_mat[imat].ppolar_matrix_tp[1][0]*spin2_fin_mvd[0]+cmc::cmc_mat[imat].ppolar_matrix_tp[1][1]*spin2_fin_mvd[1]+cmc::cmc_mat[imat].ppolar_matrix_tp[1][2]*spin2_fin_mvd[2];
			spin2_final[2]=cmc::cmc_mat[imat].ppolar_matrix_tp[2][0]*spin2_fin_mvd[0]+cmc::cmc_mat[imat].ppolar_matrix_tp[2][1]*spin2_fin_mvd[1]+cmc::cmc_mat[imat].ppolar_matrix_tp[2][2]*spin2_fin_mvd[2];

			//Calculate Energy Difference 2
			// Calculate current energy
			Eold = sim::calculate_spin_energy(atom_number2, AtomExchangeType);

         // Copy new spin position (provisionally accept move)
			atoms::x_spin_array[atom_number2] = spin2_final[0];
			atoms::y_spin_array[atom_number2] = spin2_final[1];
			atoms::z_spin_array[atom_number2] = spin2_final[2];

			// Calculate new energy
			Enew = sim::calculate_spin_energy(atom_number2, AtomExchangeType);

         // Calculate difference in Joules/mu_B
			delta_energy2 = (Enew-Eold)*mp::material[imat2].mu_s_SI*1.07828231e23; //1/9.27400915e-24

			// Calculate Delta E for both spins
			delta_energy21 = delta_energy1 + delta_energy2;

			// Compute Mz_other, Mz, Mz'
			Mz_old = cmc::cmc_mat[imat].M_other[0]*cmc::cmc_mat[imat].ppolar_vector[0] + cmc::cmc_mat[imat].M_other[1]*cmc::cmc_mat[imat].ppolar_vector[1] + cmc::cmc_mat[imat].M_other[2]*cmc::cmc_mat[imat].ppolar_vector[2];

			Mz_new = (cmc::cmc_mat[imat].M_other[0] + spin1_final[0] + spin2_final[0]- spin1_initial[0] - spin2_initial[0])*cmc::cmc_mat[imat].ppolar_vector[0] +
						(cmc::cmc_mat[imat].M_other[1] + spin1_final[1] + spin2_final[1]- spin1_initial[1] - spin2_initial[1])*cmc::cmc_mat[imat].ppolar_vector[1] +
						(cmc::cmc_mat[imat].M_other[2] + spin1_final[2] + spin2_final[2]- spin1_initial[2] - spin2_initial[2])*cmc::cmc_mat[imat].ppolar_vector[2];

			// Check for lower energy state and accept unconditionally
			//if((delta_energy21<0.0) && (Mz_new>=0.0) ) continue;
			
			// Otherwise evaluate probability for move
			//else{
				// If move is favorable then accept
				probability = exp(-delta_energy21*kBTBohr)*((Mz_new/Mz_old)*(Mz_new/Mz_old))*std::fabs(spin2_init_mvd[2]/spin2_fin_mvd[2]);
				if((probability>=mtrandom::grnd()) && (Mz_new>=0.0) ){
					cmc::cmc_mat[imat].M_other[0] = cmc::cmc_mat[imat].M_other[0] + spin1_final[0] + spin2_final[0] - spin1_initial[0] - spin2_initial[0];
					cmc::cmc_mat[imat].M_other[1] = cmc::cmc_mat[imat].M_other[1] + spin1_final[1] + spin2_final[1] - spin1_initial[1] - spin2_initial[1];
					cmc::cmc_mat[imat].M_other[2] = cmc::cmc_mat[imat].M_other[2] + spin1_final[2] + spin2_final[2] - spin1_initial[2] - spin2_initial[2];
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
				}
			//}
		}
		// if s2 not on unit sphere
		else{ 
			atoms::x_spin_array[atom_number1] = spin1_initial[0];
			atoms::y_spin_array[atom_number1] = spin1_initial[1];
			atoms::z_spin_array[atom_number1] = spin1_initial[2];
			cmc::sphere_reject+=1.0;
		}
		} // end of cmc move
		cmc::mc_total += 1.0;
	} // end of mc loop
	
	return EXIT_SUCCESS;
}

} // End of namespace sim



