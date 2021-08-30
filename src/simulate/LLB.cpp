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
///====================================================================================================
///
///       				                    	LLG
///
///  			 Subroutine to simulate an atomistic system with LLG integration scheme
///
///									Version 1.0 R Evans 02/10/2008
///
///====================================================================================================
/// \file LLG.cpp
/// Contains LLG namespace and serial version of the integrator
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "vmpi.hpp"
#include "random.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

int LLB_serial_heun(const int);

double chi_perpendicular(double x, double TC){
  //Fit Parameter
	double a0 = 0.00211549427182711;
	double a1 = 0.110224660661792;
	double a2 = -0.855153260915204;
	double a3 = 3.42088365387997;
	double a4 = -7.84821585896818;
	double a5 = 10.3247035469514;
	double a6 = -6.85608273303224;
	double a7 = 0.797453198330591;
	double a8 = 1.53787854178089;
	double a9 = -0.627128148404525;

	double chi_CGS = 0.0;
	//double chi_SI  = 0.0;
	double chi     = 0.0;
	double PI= 3.14159;

  if(x<(1.065*TC)) chi_CGS = a0+ a1*pow(pow(((1.068*TC-x)/(1.068*TC)),0.5),2.)+ a2*pow((((1.068*TC)-x)/(1.068*TC)),2.)+ a3*pow((((1.068*TC)-x)/(1.068*TC)),3.)+ a4*pow((((1.068*TC)-x)/(1.068*TC)),4.) + a5*pow((((1.068*TC)-x)/(1.068*TC)),5.) + a6*pow((((1.068*TC)-x)/(1.068*TC)),6.) + a7*pow((((1.068*TC)-x)/(1.068*TC)),7.)+ a8*pow((((1.068*TC)-x)/(1.068*TC)),8.) + a9*pow((((1.068*TC)-x)/(1.068*TC)),9.);
  else chi_CGS = (0.8*1.4/660.*TC)/(4*PI)/(x-TC);

  //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
  //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
  chi = chi_CGS*9.54393845712027; // (Tesla)

  return(chi); // [T]
}

double chi_parallel(double x, double TC){
  //Fit Parameter
  double a0 = 0.8;
  double a1 =-2.2e-07;
  double a2 = 1.95e-13;
  double a3 =-1.3e-17;
  double a4 =-4e-23;
  double a5 =-6.5076312364e-32;

  double chi_CGS = 0.0;
  //double chi_SI  = 0.0;
  double chi = 0.0;
  double PI= 3.14159;
  //double factor =  0.75947907;

  if(x<TC) chi_CGS =(a0/660.*TC)/(4.*PI)/(TC-x)+a1*pow((TC-x),1.)+ a2*pow((TC-x),3.)+a3*pow((TC-x),4.)+ a4*pow((TC-x),6.)+ a5*pow((TC-x),9.);
  else chi_CGS = (1.1*1.4/660.*TC)/(4*PI)/(x-TC);
  //chi_SI = 4*PI*chi_CGS;     // CHI_SI = 4*PI Chi_CGS
  //chi = chi_SI*4*A_FePt*A_FePt*A_FePt/MU_S_FePt/MU_0;
  chi = chi_CGS*9.54393845712027+0.308e-14; // (Tesla)

  return(chi); // [T]
}
namespace sim{
/// Master LLB Function - dispatches code path to desired LLB routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLB(const int num_steps){

   //----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){std::cout << "LLB has been called" << std::endl;}

	#ifdef MPICF
		//LLB_mpi(num_steps);
	#else
		LLB_serial_heun(num_steps);
	#endif

	return 0;
}
}



/// Performs serial Heun integration of the Landau-Lifshitz-Bloch Equation of motion
int LLB_serial_heun(const int num_steps){

	using namespace LLG_arrays;
	// Constant (to be read in later)
	const double alpha = 1.0;
	const double Tc = 661.1;
	const double temperature = sim::temperature;
	const double mu_s = 1.5E-24;
	const double n_spins = 10000.0;
	const double kB = 1.3806503e-23;
	const double gamma = 1.76e11;
	const double dt_SI = 1.0E-15;
	const double dt = dt_SI*gamma;
	const double Hx=0.0;
	const double Hy=0.0;
	const double Hz=0.0;


	// Local variables for system integration
	//const int num_atoms = atoms::num_atoms;
	double xyz[3];		// Local Delta Spin Components
	//double S_new[3];	// New Local Spin Moment
	//double mod_S;		// magnitude of spin moment

	// Setup temperature dependent variables
	double reduced_temperature = temperature/Tc;
	double Tc_o_Tc_m_T = Tc/(temperature-Tc);

	double m_e, alpha_para,alpha_perp, sigma_para, sigma_perp;
	if (temperature<=Tc){
		m_e = pow((Tc-temperature)/(Tc),0.365);
		alpha_para = alpha*(2.0/3.0)*reduced_temperature;
		alpha_perp = alpha*(1.0-temperature/(3.0*Tc));}
	else{
		m_e = 0.0;
		alpha_para = alpha*(2.0/3.0)*reduced_temperature;
		alpha_perp = alpha_para;
	}

	double m_e_squared = m_e*m_e;
	double chi_para = chi_parallel(temperature, Tc);
	double chi_perp = chi_perpendicular(temperature, Tc);;

	double one_o_2_chi_para = 1.0/(2.0*chi_para);
	double one_o_chi_perp = 1.0/chi_perp;

	if(temperature<0.1){
		sigma_para = 1.0; }
	else {
		sigma_para = sqrt(2.0*kB*temperature/(mu_s*n_spins*gamma*alpha_para*dt_SI));
	}
	sigma_perp = sqrt(2.0*kB*temperature/(mu_s*n_spins*gamma*alpha_perp*dt_SI));

	std::vector <double> Htx_perp(atoms::x_spin_array.size());
	std::vector <double> Hty_perp(atoms::x_spin_array.size());
	std::vector <double> Htz_perp(atoms::x_spin_array.size());
	std::vector <double> Htx_para(atoms::x_spin_array.size());
	std::vector <double> Hty_para(atoms::x_spin_array.size());
	std::vector <double> Htz_para(atoms::x_spin_array.size());

	// precalculate thermal fields
	generate (Htx_perp.begin(),Htx_perp.end(), mtrandom::gaussian);
	generate (Hty_perp.begin(),Hty_perp.end(), mtrandom::gaussian);
	generate (Htz_perp.begin(),Htz_perp.end(), mtrandom::gaussian);
	generate (Htx_para.begin(),Htx_para.end(), mtrandom::gaussian);
	generate (Hty_para.begin(),Hty_para.end(), mtrandom::gaussian);
	generate (Htz_para.begin(),Htz_para.end(), mtrandom::gaussian);

	for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
		Htx_perp[atom] *= sigma_perp;
		Hty_perp[atom] *= sigma_perp;
		Htz_perp[atom] *= sigma_perp;
		Htx_para[atom] *= sigma_para;
		Hty_para[atom] *= sigma_para;
		Htz_para[atom] *= sigma_para;
	}

	for(int t=0;t<num_steps;t++){

		// Store initial spin positions
		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
			x_initial_spin_array[atom] = atoms::x_spin_array[atom];
			y_initial_spin_array[atom] = atoms::y_spin_array[atom];
			z_initial_spin_array[atom] = atoms::z_spin_array[atom];
		}

		// Set field arrays to zero
		fill (atoms::x_total_spin_field_array.begin(),atoms::x_total_spin_field_array.end(),0.0);
		fill (atoms::y_total_spin_field_array.begin(),atoms::y_total_spin_field_array.end(),0.0);
		fill (atoms::z_total_spin_field_array.begin(),atoms::z_total_spin_field_array.end(),0.0);
		fill (atoms::x_total_external_field_array.begin(),atoms::x_total_external_field_array.end(),Hx);
		fill (atoms::y_total_external_field_array.begin(),atoms::y_total_external_field_array.end(),Hy);
		fill (atoms::z_total_external_field_array.begin(),atoms::z_total_external_field_array.end(),Hz);

		// Calculate fields
		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
			double m[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];
			double pf;
			if(temperature<=Tc){
				pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
			}
			else{
				pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);
			}

			atoms::x_total_spin_field_array[atom] = (pf-one_o_chi_perp)*m[0];
			atoms::y_total_spin_field_array[atom] = (pf-one_o_chi_perp)*m[1];
			atoms::z_total_spin_field_array[atom] = (pf-0.0				 )*m[2];
		}

		// Calculate Euler Step
		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};
			const double H_perp[3]={H[0]+Htx_perp[atom], H[1]+Hty_perp[atom], H[2]+Htz_perp[atom]};
			const double H_para[3]={H[0]+Htx_para[atom], H[1]+Hty_para[atom], H[2]+Htz_para[atom]};
			const double one_o_m_squared = 1.0/(S[1]*S[1]+S[2]*S[2]+S[0]*S[0]);

			// Calculate Delta S
			xyz[0]= 	-(S[1]*H[2]-S[2]*H[1])
						+ alpha_para*S[0]*S[0]*H_para[0]*one_o_m_squared
						-alpha_perp*(S[1]*(S[0]*H_perp[1]-S[1]*H_perp[0])-S[2]*(S[2]*H_perp[0]-S[0]*H_perp[2]))*one_o_m_squared;

			xyz[1]= 	-(S[2]*H[0]-S[0]*H[2])
						+ alpha_para*S[1]*S[1]*H_para[1]*one_o_m_squared
						-alpha_perp*(S[2]*(S[1]*H_perp[2]-S[2]*H_perp[1])-S[0]*(S[0]*H_perp[1]-S[1]*H_perp[0]))*one_o_m_squared;

			xyz[2]=	-(S[0]*H[1]-S[1]*H[0])
						+ alpha_para*S[2]*S[2]*H_para[2]*one_o_m_squared
						-alpha_perp*(S[0]*(S[2]*H_perp[0]-S[0]*H_perp[2])-S[1]*(S[1]*H_perp[2]-S[2]*H_perp[1]))*one_o_m_squared;

			// Store dS in euler array
			x_euler_array[atom]=xyz[0];
			y_euler_array[atom]=xyz[1];
			z_euler_array[atom]=xyz[2];

			// Calculate Euler Step
			x_spin_storage_array[atom]=S[0]+xyz[0]*dt;
			y_spin_storage_array[atom]=S[1]+xyz[1]*dt;
			z_spin_storage_array[atom]=S[2]+xyz[2]*dt;

		}

		// Copy new spins to spin array
		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
			atoms::x_spin_array[atom]=x_spin_storage_array[atom];
			atoms::y_spin_array[atom]=y_spin_storage_array[atom];
			atoms::z_spin_array[atom]=z_spin_storage_array[atom];
		}

		// Recalculate spin dependent fields
		fill (atoms::x_total_spin_field_array.begin(),atoms::x_total_spin_field_array.end(),0.0);
		fill (atoms::y_total_spin_field_array.begin(),atoms::y_total_spin_field_array.end(),0.0);
		fill (atoms::z_total_spin_field_array.begin(),atoms::z_total_spin_field_array.end(),0.0);

		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
			double m[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			double m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];
			double pf;
			if(temperature<=Tc){
				pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
			}
			else{
				pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);
			}

			atoms::x_total_spin_field_array[atom] = (pf-one_o_chi_perp)*m[0];
			atoms::y_total_spin_field_array[atom] = (pf-one_o_chi_perp)*m[1];
			atoms::z_total_spin_field_array[atom] = (pf-0.0				 )*m[2];
		}

		// Calculate Heun Gradients

		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){

			// Store local spin in Sand local field in H
			const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
			const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
										atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
										atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};
			const double H_perp[3]={H[0]+Htx_perp[atom], H[1]+Hty_perp[atom], H[2]+Htz_perp[atom]};
			const double H_para[3]={H[0]+Htx_para[atom], H[1]+Hty_para[atom], H[2]+Htz_para[atom]};
			const double one_o_m_squared = 1.0/(S[1]*S[1]+S[2]*S[2]+S[0]*S[0]);

			// Calculate Delta S
			xyz[0]= 	-(S[1]*H[2]-S[2]*H[1])
						+ alpha_para*S[0]*S[0]*H_para[0]*one_o_m_squared
						-alpha_perp*(S[1]*(S[0]*H_perp[1]-S[1]*H_perp[0])-S[2]*(S[2]*H_perp[0]-S[0]*H_perp[2]))*one_o_m_squared;

			xyz[1]= 	-(S[2]*H[0]-S[0]*H[2])
						+ alpha_para*S[1]*S[1]*H_para[1]*one_o_m_squared
						-alpha_perp*(S[2]*(S[1]*H_perp[2]-S[2]*H_perp[1])-S[0]*(S[0]*H_perp[1]-S[1]*H_perp[0]))*one_o_m_squared;

			xyz[2]=	-(S[0]*H[1]-S[1]*H[0])
						+ alpha_para*S[2]*S[2]*H_para[2]*one_o_m_squared
						-alpha_perp*(S[0]*(S[2]*H_perp[0]-S[0]*H_perp[2])-S[1]*(S[1]*H_perp[2]-S[2]*H_perp[1]))*one_o_m_squared;

			// Store dS in heun array
			x_heun_array[atom]=xyz[0];
			y_heun_array[atom]=xyz[1];
			z_heun_array[atom]=xyz[2];
		}

		//----------------------------------------
		// Calculate Heun Step
		//----------------------------------------

		for(unsigned int atom=0;atom<atoms::x_spin_array.size();atom++){
			atoms::x_spin_array[atom]=x_initial_spin_array[atom]+0.5*dt*(x_euler_array[atom]+x_heun_array[atom]);
			atoms::y_spin_array[atom]=y_initial_spin_array[atom]+0.5*dt*(y_euler_array[atom]+y_heun_array[atom]);
			atoms::z_spin_array[atom]=z_initial_spin_array[atom]+0.5*dt*(z_euler_array[atom]+z_heun_array[atom]);
		}

	}

	return EXIT_SUCCESS;
	}
