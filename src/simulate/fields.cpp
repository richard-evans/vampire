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
//=============================================================================
//
//                                     Fields
//
//              Subroutines to calculate fields for the hamiltonian
//
//                         Version 1.0 R Evans 20/10/2008
//
//====================================================================================================
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "environment.hpp"
#include "dipole.hpp"
#include "hamr.hpp"
#include "ltmp.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "spintorque.hpp"
#include "spintransport.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "../micromagnetic/internal.hpp"

// sim module header
#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

//========================
//function prototypes
//========================

int calculate_exchange_fields(const int,const int);
int calculate_applied_fields(const int,const int);
int calculate_thermal_fields(const int,const int);
int calculate_dipolar_fields(const int,const int);
void calculate_fmr_fields(const int,const int);
void calculate_lagrange_fields(const int,const int);
void calculate_full_spin_fields(const int start_index,const int end_index);

namespace sim{

void calculate_spin_fields(const int start_index,const int end_index){

	///======================================================
	/// 		Subroutine to calculate spin dependent fields
	///
	///			Version 1.0 R Evans 20/10/2008
	///======================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_spin_fields has been called" << std::endl;}

	// Initialise Total Spin Fields to zero
	fill (atoms::x_total_spin_field_array.begin()+start_index,atoms::x_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_spin_field_array.begin()+start_index,atoms::y_total_spin_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_spin_field_array.begin()+start_index,atoms::z_total_spin_field_array.begin()+end_index,0.0);

   //-----------------------------------------
	// Calculate exchange Fields
   //-----------------------------------------
   exchange::fields(start_index, // first atom for exchange interactions to be calculated
                    end_index, // last +1 atom to be calculated
                    atoms::neighbour_list_start_index,
                    atoms::neighbour_list_end_index,
                    atoms::type_array, // type for atom
                    atoms::neighbour_list_array, // list of interactions between atoms
                    atoms::neighbour_interaction_type_array, // list of interaction type for each pair of atoms with value given in exchange list
                    atoms::i_exchange_list, // list of isotropic exchange constants
                    atoms::v_exchange_list, // list of vectorial exchange constants
                    atoms::t_exchange_list, // list of tensorial exchange constants
                    atoms::x_spin_array,
                    atoms::y_spin_array,
                    atoms::z_spin_array,
                    atoms::x_total_spin_field_array,
                    atoms::y_total_spin_field_array,
                    atoms::z_total_spin_field_array);

   //-----------------------------------------
   // calculate anistropy fields
   //-----------------------------------------
   anisotropy::fields(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::type_array,
                      atoms::x_total_spin_field_array, atoms::y_total_spin_field_array, atoms::z_total_spin_field_array,
                      start_index, end_index, sim::temperature);

	// Spin Dependent Extra Fields
	if(sim::lagrange_multiplier==true) calculate_lagrange_fields(start_index,end_index);

	// Add spin torque fields
	if( sim::internal::enable_spin_torque_fields == true ||
		 sim::internal::enable_vcma_fields        == true ){
		calculate_full_spin_fields(start_index,end_index);
	}

	return;

}


void calculate_external_fields(const int start_index,const int end_index){
	///======================================================
	/// 		Subroutine to calculate external fields
	///
	///			Version 1.0 R Evans 20/10/2008
	///======================================================
	//const int num_atoms = atoms::num_atoms;

	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "calculate_external_fields has been called" << std::endl;}

	// Initialise Total External Fields to zero
	fill (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index,0.0);
	fill (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index,0.0);

	if(program::program==7){

		// Calculate thermal field and applied field due to HAMR process
		hamr::fields(start_index,
					end_index,
					sim::H_applied,
					sim::temperature,
					sim::H_vec[0],
					sim::H_vec[1],
					sim::H_vec[2],
					atoms::x_total_external_field_array,
					atoms::y_total_external_field_array,
					atoms::z_total_external_field_array);
	}
   else if(program::program==13){

      // Local thermal Fields
      ltmp::get_localised_thermal_fields(atoms::x_total_external_field_array,atoms::y_total_external_field_array,
            atoms::z_total_external_field_array, start_index, end_index);

      // Applied Fields
      if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);

   }
	else{

		// Thermal Fields
		if(sim::hamiltonian_simulation_flags[3]==1) calculate_thermal_fields(start_index,end_index);

		// Applied Fields
		if(sim::hamiltonian_simulation_flags[2]==1) calculate_applied_fields(start_index,end_index);

	}

   // Get updated spin torque fields
   st::get_spin_torque_fields(atoms::x_total_external_field_array, atoms::y_total_external_field_array, atoms::z_total_external_field_array, start_index, end_index);

   // Get updated spin torque fields
   spin_transport::calculate_field(start_index, end_index, atoms::x_total_external_field_array, atoms::y_total_external_field_array, atoms::z_total_external_field_array);

	// FMR Fields only for fmr program
	if(sim::enable_fmr) calculate_fmr_fields(start_index,end_index);

	// Dipolar Fields
	calculate_dipolar_fields(start_index,end_index);

	return;
}
}

int calculate_applied_fields(const int start_index,const int end_index){
	///==========================================================================
	///
	/// 	Function to calculate applied fields
	///
	///		Version 1.0 R Evans 20/10/2008
	///		Version 2.0 R F L Evans 18/11/2012
	///
	///==========================================================================

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "calculate_applied_fields has been called" << std::endl;}

	// Declare constant temporaries for global field
	const double Hx=sim::H_vec[0]*sim::H_applied;
	const double Hy=sim::H_vec[1]*sim::H_applied;
	const double Hz=sim::H_vec[2]*sim::H_applied;

	// Declare array for local (material specific) applied field
	std::vector<double> Hlocal(0);

	// Check for local applied field
	if(sim::local_applied_field==true){
		Hlocal.reserve(3*mp::material.size());

		// Loop over all materials
		for(unsigned int mat=0;mat<mp::material.size();mat++){
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[0]);
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[1]);
			Hlocal.push_back(mp::material[mat].applied_field_strength*mp::material[mat].applied_field_unit_vector[2]);
		}

		// Add local field AND global field
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			atoms::x_total_external_field_array[atom] += Hx + Hlocal[3*imaterial + 0];
			atoms::y_total_external_field_array[atom] += Hy + Hlocal[3*imaterial + 1];
			atoms::z_total_external_field_array[atom] += Hz + Hlocal[3*imaterial + 2];
		}
	}
	else{
		// Calculate global field
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}

	}

	// Add external field from thin film sample
	if(sim::ext_demag==true){

      const std::vector<double> m_l = stats::system_magnetization.get_magnetization();

		// calculate global demag field -mu_0 M D, M = m/V
		const double mu_0= -4.0*M_PI*1.0e-7/(cs::system_dimensions[0]*cs::system_dimensions[1]*cs::system_dimensions[2]*1.0e-30);
      const double HD[3]={	mu_0*sim::demag_factor[0]*m_l[0],
                           mu_0*sim::demag_factor[1]*m_l[1],
                           mu_0*sim::demag_factor[2]*m_l[2]};

		//std::cout << "mu_0" << "\t" << mu_0 << std::endl;
		//std::cout << "Magnetisation " << stats::total_mag_actual[0] << "\t" << stats::total_mag_actual[1] << "\t" << stats::total_mag_actual[2] << std::endl;
		//std::cout << "External Demag Field " << HD[0] << "\t" << HD[1] << "\t" << HD[2] << std::endl;
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += HD[0];
			atoms::y_total_external_field_array[atom] += HD[1];
			atoms::z_total_external_field_array[atom] += HD[2];
		}
	}

	if(micromagnetic::internal::bias_magnets == true){
				for(int atom=start_index;atom<end_index;atom++){
					atoms::x_total_external_field_array[atom] += micromagnetic::atomistic_bias_field_x[atom];
					atoms::y_total_external_field_array[atom] += micromagnetic::atomistic_bias_field_y[atom];
					atoms::z_total_external_field_array[atom] += micromagnetic::atomistic_bias_field_z[atom];
			//		std::cout << atom << '\t' << micromagnetic::atomistic_bias_field_x[atom] << '\t' << std::endl;
				}
	}
	if(environment::enabled == true){
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += environment::atomistic_environment_field_x[atom];
			atoms::y_total_external_field_array[atom] += environment::atomistic_environment_field_y[atom];
			atoms::z_total_external_field_array[atom] += environment::atomistic_environment_field_z[atom];
	//		std::cout << atom << '\t' << micromagnetic::atomistic_bias_field_x[atom] << '\t' << std::endl;
		}
	}


	return 0;

}

int calculate_thermal_fields(const int start_index,const int end_index){
   ///======================================================
   /// 		Subroutine to calculate thermal fields
   ///
   ///      Version 1.2 R Evans 12/08/2014
   ///======================================================

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "calculate_thermal_fields has been called" << std::endl;}

   // unroll sigma for speed
   std::vector<double> sigma_prefactor(0);
   sigma_prefactor.reserve(mp::material.size());

   // Calculate material temperature (with optional rescaling)
   for(unsigned int mat=0;mat<mp::material.size();mat++){
      double temperature = sim::temperature;
      // Check for localised temperature
      if(sim::local_temperature) temperature = mp::material[mat].temperature;
      // Calculate temperature rescaling
      double alpha = mp::material[mat].temperature_rescaling_alpha;
      double Tc = mp::material[mat].temperature_rescaling_Tc;
      // if T<Tc T/Tc = (T/Tc)^alpha else T = T
      double rescaled_temperature = temperature < Tc ? Tc*pow(temperature/Tc,alpha) : temperature;
      double sqrt_T=sqrt(rescaled_temperature);
      sigma_prefactor.push_back(sqrt_T*mp::material[mat].H_th_sigma);
   }

   generate (atoms::x_total_external_field_array.begin()+start_index,atoms::x_total_external_field_array.begin()+end_index, mtrandom::gaussian);
   generate (atoms::y_total_external_field_array.begin()+start_index,atoms::y_total_external_field_array.begin()+end_index, mtrandom::gaussian);
   generate (atoms::z_total_external_field_array.begin()+start_index,atoms::z_total_external_field_array.begin()+end_index, mtrandom::gaussian);

   for(int atom=start_index;atom<end_index;atom++){

      const int imaterial=atoms::type_array[atom];
      const double H_th_sigma = sigma_prefactor[imaterial];

      atoms::x_total_external_field_array[atom] *= H_th_sigma;
		atoms::y_total_external_field_array[atom] *= H_th_sigma;
		atoms::z_total_external_field_array[atom] *= H_th_sigma;
	}

   return EXIT_SUCCESS;
}

int calculate_dipolar_fields(const int start_index,const int end_index){

	///======================================================
	/// 		Subroutine to calculate dipolar fields
	///
	///			Version 1.0 R Evans 02/11/2009
	///======================================================
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
   if(err::check==true){std::cout << "calculate_dipolar_fields has been called" << std::endl;}

   // Add dipolar fields
   if(dipole::activated){
      for(int atom=start_index;atom<end_index;atom++){
         atoms::x_total_external_field_array[atom] += dipole::atom_dipolar_field_array_x[atom];
         atoms::y_total_external_field_array[atom] += dipole::atom_dipolar_field_array_y[atom];
         atoms::z_total_external_field_array[atom] += dipole::atom_dipolar_field_array_z[atom];
         /*std::cout << atoms::x_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_x[atom] << "\t";
         std::cout << atoms::y_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_y[atom] << "\t";
         std::cout << atoms::z_total_external_field_array[atom] << "\t" <<  dipole::atom_dipolar_field_array_z[atom] << std::endl;*/
      }
   }

   return 0;
}

void calculate_fmr_fields(const int start_index,const int end_index){

	if(err::check==true){std::cout << "calculate_fmr_fields has been called" << std::endl;}

	// Calculate fmr constants
	const double real_time = sim::time*mp::dt_SI;
	const double omega = sim::fmr_field_frequency; // Hz
	const double Hfmrx = sim::fmr_field_unit_vector[0];
	const double Hfmry = sim::fmr_field_unit_vector[1];
	const double Hfmrz = sim::fmr_field_unit_vector[2];
	const double Hsinwt = sim::fmr_field_strength * sin(2.0 * M_PI * omega * real_time);
	const double Hx = Hfmrx * Hsinwt;
	const double Hy = Hfmry * Hsinwt;
	const double Hz = Hfmrz * Hsinwt;

	// Save fmr field strength for possible output
	sim::fmr_field = Hsinwt;

	if(sim::local_fmr_field==true){

		std::vector<double> H_fmr_local;
		H_fmr_local.reserve(3*mp::material.size());

		// Loop over all materials
		for(unsigned int mat=0;mat<mp::material.size();mat++){
			const double Hsinwt_local = mp::material[mat].fmr_field_strength * sin( 2.0 * M_PI * real_time * mp::material[mat].fmr_field_frequency );

			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[0]);
			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[1]);
			H_fmr_local.push_back(Hsinwt_local*mp::material[mat].fmr_field_unit_vector[2]);

		}

		// Add local field AND global field
		for(int atom=start_index;atom<end_index;atom++){
			const int imaterial=atoms::type_array[atom];
			atoms::x_total_external_field_array[atom] += Hx + H_fmr_local[3*imaterial + 0];
			atoms::y_total_external_field_array[atom] += Hy + H_fmr_local[3*imaterial + 1];
			atoms::z_total_external_field_array[atom] += Hz + H_fmr_local[3*imaterial + 2];
		}
	}
	else{
		// Add fmr field
		for(int atom=start_index;atom<end_index;atom++){
			atoms::x_total_external_field_array[atom] += Hx;
			atoms::y_total_external_field_array[atom] += Hy;
			atoms::z_total_external_field_array[atom] += Hz;
		}
	}

	return;

}

///------------------------------------------------------
///  Function to calculate LaGrange multiplier fields for
///  constrained minimization
///
///  (c) R F L Evans 2013
///
///------------------------------------------------------
void calculate_lagrange_fields(const int start_index,const int end_index){

   // LaGrange Multiplier
   const double lx=sim::lagrange_lambda_x;
   const double ly=sim::lagrange_lambda_y;
   const double lz=sim::lagrange_lambda_z;

   // Constraint vector
   const double nu_x=cos(sim::constraint_theta*M_PI/180.0)*sin(sim::constraint_phi*M_PI/180.0);
   const double nu_y=sin(sim::constraint_theta*M_PI/180.0)*sin(sim::constraint_phi*M_PI/180.0);
   const double nu_z=cos(sim::constraint_phi*M_PI/180.0);

   // Magnetisation
   const double imm=1.0/sim::lagrange_m;
   const double imm3=1.0/(sim::lagrange_m*sim::lagrange_m*sim::lagrange_m);

   const double N=sim::lagrange_N;

   // Calculate LaGrange fields
   for(int atom=start_index;atom<end_index;atom++){
      const double sx=atoms::x_spin_array[atom];
      const double sy=atoms::y_spin_array[atom];
      const double sz=atoms::z_spin_array[atom];

      //std::cout << "S " << sx << "\t" << sy << "\t" << sz << std::endl;
      //std::cout << "L " << lx << "\t" << ly << "\t" << lz << std::endl;
      //std::cout << imm << "\t" << imm3 << std::endl;

      const double lambda_dot_s = lx*sx + ly*sy + lz*sz;

      atoms::x_total_spin_field_array[atom]+=N*(lx*imm - lambda_dot_s*sx*imm3 - nu_x);
      atoms::y_total_spin_field_array[atom]+=N*(ly*imm - lambda_dot_s*sy*imm3 - nu_y);
      atoms::z_total_spin_field_array[atom]+=N*(lz*imm - lambda_dot_s*sz*imm3 - nu_z);

      //std::cout << "\t" << N*(lx*imm - lambda_dot_s*sx*imm3 - nu_x) << std::endl;
      //std::cout << "\t" << N*(ly*imm - lambda_dot_s*sy*imm3 - nu_y) << std::endl;
      //std::cout << "\t" << N*(lz*imm - lambda_dot_s*sz*imm3 - nu_z) << std::endl;
      //std::cin.get();
   }
   return;

}

//------------------------------------------------------------------------------
// Master function to calculate fields in large loop
//------------------------------------------------------------------------------
void calculate_full_spin_fields(const int start_index,const int end_index){

	using namespace sim::internal;

   for(int atom=start_index;atom<end_index;atom++){

		// temporary variables for field components
		double hx = 0.0;
		double hy = 0.0;
		double hz = 0.0;

		// temporary constant for spin components
		const double sx = atoms::x_spin_array[atom];
		const double sy = atoms::y_spin_array[atom];
		const double sz = atoms::z_spin_array[atom];

		// get material parameter
		const int material = atoms::type_array[atom];

		const double alpha = mp::material[material].alpha;
		//----------------------------------------------------------------------------------
		// Slonczewski spin torque field
		//----------------------------------------------------------------------------------

		// save polarization to temporary constant
		const double stpx = stt_polarization_unit_vector[0];
		const double stpy = stt_polarization_unit_vector[1];
		const double stpz = stt_polarization_unit_vector[2];

		const double strj = stt_rj[material];
		const double stpj = stt_pj[material];

		const double stt_lambda = stt_asm[material];
		const double factor = program::fractional_electric_field_strength / (1.0 + stt_lambda*(sx*stpx + sy*stpy + sz*stpz) );

		// calculate field
		hx += factor * ( (strj-alpha*stpj)*(sy*stpz - sz*stpy) + (stpj+alpha*strj)*stpx );
		hy += factor * ( (strj-alpha*stpj)*(sz*stpx - sx*stpz) + (stpj+alpha*strj)*stpy );
		hz += factor * ( (strj-alpha*stpj)*(sx*stpy - sy*stpx) + (stpj+alpha*strj)*stpz );

		//----------------------------------------------------------------------------------
		// Spin orbit torque (SOT) field
		//----------------------------------------------------------------------------------

		// save polarization to temporary constant
		const double sotpx = sot_polarization_unit_vector[0];
		const double sotpy = sot_polarization_unit_vector[1];
		const double sotpz = sot_polarization_unit_vector[2];

		const double sotrj = sot_rj[material];
		const double sotpj = sot_pj[material];

		const double sot_lambda = sot_asm[material];
		const double sot_factor = program::fractional_electric_field_strength / (1.0 + sot_lambda*(sx*sotpx + sy*sotpy + sz*sotpz) );

		// calculate field
		hx += sot_factor * ( (sotrj-alpha*sotpj)*(sy*sotpz - sz*sotpy) + (sotpj+alpha*sotrj)*sotpx );
		hy += sot_factor * ( (sotrj-alpha*sotpj)*(sz*sotpx - sx*sotpz) + (sotpj+alpha*sotrj)*sotpy );
		hz += sot_factor * ( (sotrj-alpha*sotpj)*(sx*sotpy - sy*sotpx) + (sotpj+alpha*sotrj)*sotpz );

		//----------------------------------------------------------------------------------
		// VCMA field
		//----------------------------------------------------------------------------------
		const double vcma = program::fractional_electric_field_strength * spin_transport::get_voltage() * sim::internal::vcmak[material];
		hz += 2.0 * vcma * sz;

		//----------------------------------------------------------------------------------
		// save field to spin field array
		//----------------------------------------------------------------------------------
		atoms::x_total_spin_field_array[atom] += hx;
		atoms::y_total_spin_field_array[atom] += hy;
		atoms::z_total_spin_field_array[atom] += hz;

	}

	return;

}
