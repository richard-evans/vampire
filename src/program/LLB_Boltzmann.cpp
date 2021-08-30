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
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace program{
	
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

int LLB_Boltzmann(){
  // check calling of routine if error checking is activated
  if(err::check==true){std::cout << "program::LLB_Boltzmann has been called" << std::endl;}


	// Initialise spins to z-direction
	for(int atom =0;atom<atoms::num_atoms;atom++){
			atoms::x_spin_array[atom]=0.0;
			atoms::y_spin_array[atom]=0.0;
			atoms::z_spin_array[atom]=1.0;
	}
	
	for(int T=300;T<310;T+=10){
		sim::temperature=double(T);
		double mean_M=0.0;
		double counter=0.0;
		
		double P[101][101];
		double P1D[1001];
		for(int para=0;para<101;para++){
			for(int perp=0;perp<101;perp++){
				P[para][perp]=0.0;
			}
		}
		for(int para=0;para<1001;para++){
			P1D[para]=0.0;
		}
		// Simulate system
		for(sim::time=0;sim::time<10000000;sim::time+=1){
			// Calculate LLB
			sim::LLB(1);
			  if(sim::time%100000==0){
				  std::cout << sim::time << std::endl;
			  }
			// Calculate mag_m, mag
			if(sim::time>10000){
				double S[3] = {atoms::x_spin_array[0],atoms::y_spin_array[0],atoms::z_spin_array[0]};
				double mag_m=sqrt(S[0]*S[0]+S[1]*S[1]+S[2]*S[2]);
				double mz=S[2];
				double mx=sqrt(S[0]*S[0]); 
				int para = int(mz*100.0+0.5);
				int perp = int(mx*100.0+0.5);
				int para1D = int(mag_m*1000.0+0.5);
				P[para][perp]+=1.0;
				P1D[para1D]+=1.0;
				mean_M+=mag_m;
				counter+=1.0;
			}
		}
		std::cout << sim::temperature << "\t" << atoms::x_spin_array[0];
		std::cout << "\t" << atoms::y_spin_array[0];
		std::cout << "\t" << atoms::z_spin_array[0];
		std::cout << "\t" << mean_M/counter;
		std::cout << std::endl;
		zmag << sim::temperature << "\t" << mean_M/counter << std::endl;
		std::ofstream pfile("LLBprob");
		const double Tc = 661.1;
		const double chi_para = chi_parallel(sim::temperature, Tc);
		const double chi_perp = chi_perpendicular(sim::temperature, Tc);
		const double n_spins = 10000.0;
		const double kB = 1.3806503e-23;
		const double mu_s = 1.5E-24;
		for(int para=0;para<101;para++){
			for(int perp=0;perp<101;perp++){
			double mp=double(para)/100.0;
			double mt=double(perp)/100.0;
			double m_e = pow((Tc-sim::temperature)/(Tc),0.365);
			double F = n_spins*mu_s*(	((mp*mp-m_e*m_e)*(mp*mp-m_e*m_e))/(8.0*chi_para*m_e*m_e)	+ mt*mt/(2.0*chi_perp));
			double PF = exp(-F/(kB*sim::temperature));
				pfile << double(para)/100.0 << "\t" << double(perp)/100.0 << "\t" << P[para][perp]/counter << "\t" << PF << std::endl;
			}
			pfile << std::endl;
		}
		//===========================================================
		std::ofstream pfile1D("LLBprob1D");
		//const double Tc = 661.1;
		//const double chi_para = chi_parallel(sim::temperature, Tc);
		//const double n_spins = 10000.0;
		//const double kB = 1.3806503e-23;
		//const double mu_s = 1.5E-24;
		std::cout << "m_e: " << pow((Tc-sim::temperature)/(Tc),0.365) << std::endl;
		for(int para=0;para<1001;para++){
			double m=double(para)/1000.0;
			double m_e = pow((Tc-sim::temperature)/(Tc),0.365);
			double F = n_spins*mu_s*(	((m*m-m_e*m_e)*(m*m-m_e*m_e))/(8.0*chi_para*m_e*m_e)	);
			double PF = exp(-F/(kB*sim::temperature));
			//std::cout << "F: " << F << " PF: " << PF << std::endl;
			pfile1D << (double(para))/1000.0 << "\t" << P1D[para]/counter << "\t" << PF << std::endl;
		}

	}

	return EXIT_SUCCESS;
}


}//end of namespace program


