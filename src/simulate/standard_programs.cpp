///
/// @file
/// @brief This is the brief (one line only) description of the funtion of this file. 
///
/// @details This is the detailed description of the funtion of this file
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
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    11/01/2010
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
#include "atoms.hpp"
#include "material.hpp"
#include "public.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include <sstream>
#include <iomanip>
#include <iostream>

	int LLG(int const);
	int LLG_relax(int const);
	int initialise_system();
	
/// @namespace program
/// @brief A Namespace containing functions for standard programs.
/// 
/// @internal
///=====================================================================================
///
namespace program{

/// @brief Function to calculate the temperature dependence of the magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// accoring to the input flag - either randomly or ordered.For the ordered case the temperature sequence
/// increases from zero, for the random case the temperature decreases from the maximum temperature. After
/// initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
/// @section notes Implementation Notes
/// Capable of hot>cold or cold>hot calculation. 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    11/01/2010
///
/// @param[in] init Determines whether the system is initialised randomly (0) or ordered (1)
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		11/01/2010
///	Revision:	  ---
///=====================================================================================
///
int curie_temperature(bool init){
	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){std::cout << "program::curie_temperature has been called" << std::endl;}
	
	// function prototypes

	//int output_pov_file();

	// Setup LLG arrays
	sim::initialise();

	// Initialise spins to random state
	for(int atom =0;atom<atoms::num_atoms;atom++){
		if(init==false){
			atoms::x_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			atoms::y_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
			atoms::z_spin_array[atom]=2.0*mtrandom::grnd()-1.0;
		}
		else{
			double parity = double(atoms::grain_array[atom]%2);
			atoms::x_spin_array[atom]=0.0;			
			atoms::y_spin_array[atom]=0.0;			
			atoms::z_spin_array[atom]=2.0*parity-1.0;
		}
	}
		//sim::H_applied=0.0;
	// Set up loop variables
	
	vout::pov_file();
		
	//      Perform Temperature Loop
	for(int temperature=0;temperature<=1000;temperature+=10){
		// Set system temperature
		sim::temperature=double(temperature); 
		
		// Equilibrate system
		LLG(sim::equilibration_time);
		
		// Simulate system
		for(sim::time=0;sim::time<sim::loop_time;sim::time+=sim::partial_time){
			
			// Calculate LLG
			LLG(sim::partial_time);
			
			// Calculate mag_m, mag
			stats::mag_m();
		}
		// Output to screen and file after each temperature
		if(vmpi::my_rank==0){
			std::cout << sim::temperature << "\t" << stats::total_mag_m_norm;
			std::cout << "\t" << stats::total_mag_norm[0];
			std::cout << "\t" << stats::total_mag_norm[1];
			std::cout << "\t" << stats::total_mag_norm[2];
			std::cout << std::endl;
			vmag << sim::temperature << "\t" << stats::total_mag_m_norm << std::endl;
		}
		
	vout::pov_file();
		
	} // End of temperature loop


		
	return EXIT_SUCCESS;
}



int bmark(){
  // check calling of routine if error checking is activated
  if(error_checking::error_check==true){std::cout << "program::bmark has been called" << std::endl;}

  // Setup LLG arrays
  sim::initialise();

	// Initialise spins to random state
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
  }

  sim::temperature=300.0;

  // Simulate system
  for(sim::time=0;sim::time<10000;sim::time+=1){

  // Calculate LLG
  LLG(1);

      // Calculate mag_m, mag
  if(sim::time%1000==0){
      stats::mag_m();
		//vout::pov_file();
  if(vmpi::my_rank==0){
    std::cout << sim::time << "\t" << stats::total_mag_m_norm;
    std::cout << "\t" << stats::total_mag_norm[0];
    std::cout << "\t" << stats::total_mag_norm[1];
    std::cout << "\t" << stats::total_mag_norm[2];
    std::cout << std::endl;
    //vmag << sim::temperature << "\t" << stats::total_mag_m_norm << std::endl;
  }

  }
  
  }



return EXIT_SUCCESS;
}
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
  if(error_checking::error_check==true){std::cout << "program::LLB_Boltzmann has been called" << std::endl;}

  // Setup LLG arrays
  sim::initialise();

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
			// Calculate LLG
			sim::LLB(1);
			  if(sim::time%100000==0){
				  std::cout << sim::time << std::endl;
			  }
			// Calculate mag_m, mag
			if(sim::time>10000){
				double S[3] = {atoms::x_spin_array[0],atoms::y_spin_array[0],atoms::z_spin_array[0]};
				double mag_m=sqrt(S[0]*S[0]+S[1]*S[1]+S[2]*S[2]);
				double mz=S[2];
				double mx=sqrt(S[0]*S[0]); //sqrt(S[0]*S[0]+S[1]*S[1]);
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
		vmag << sim::temperature << "\t" << mean_M/counter << std::endl;
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

/// @brief Function to calculate a time series with gaussian cooling profile
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a time sequence of sub-calculations of fixed temperature. The system is initialised 
/// ordered. After initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    27/01/2010
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		27/01/2010
///	Revision:	  ---
///=====================================================================================
///
int hamr_run(){
	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){std::cout << "program::hamr_run has been called" << std::endl;}
	
	// function prototypes


	// Setup LLG arrays
	sim::initialise();

	// Initialise spins to +z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
	}
	
	// Set up loop variables
	sim::H_applied=-0.8;
	
	double cooling_time=2.0e-9; // Seconds
	double max_dT = 300.0; 
	double RT = 300.0;
	// Set initial system temperature
	sim::temperature=300.0;
	
	// Equilibrate system
	LLG(sim::equilibration_time);
	
	// Simulate system with single timestep resolution
	for(sim::time=0;sim::time<sim::loop_time;sim::time++){
		
		// calculate real time and temperature using gaussian cooling
		double actual_time = double(sim::time)*mp::dt_SI;
		sim::temperature=RT+max_dT*exp(-((actual_time*actual_time)/(cooling_time*cooling_time)));
		
		// Calculate LLG
		LLG(1);
		
		// Calculate mag_m, mag
		if(sim::time%sim::partial_time==0){
			stats::mag_m();
			// Output to screen and file after each temperature
			if(vmpi::my_rank==0){
			  //std::cout <<
			  //std::cout << sim::temperature << "\t" << stats::total_mag_m_norm;
			  //std::cout << "\t" << stats::total_mag_norm[0];
			  //std::cout << "\t" << stats::total_mag_norm[1];
			  //std::cout << "\t" << stats::total_mag_norm[2];
			  //std::cout << std::endl;
			  vmag << actual_time << "\t" << sim::temperature << "\t" << stats::total_mag_m_norm;
				vmag << "\t" << stats::total_mag_norm[0];
				vmag << "\t" << stats::total_mag_norm[1];
				vmag << "\t" << stats::total_mag_norm[2];
				vmag << std::endl;
			}
		}
	}
	//vout::pov_file();
		
	return EXIT_SUCCESS;
}

/// @brief Function to calculate a time series with two temperature model heating/cooling profile
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a time sequence of sub-calculations of fixed temperature. The system is initialised 
/// ordered. After initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    10/03/2010
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		10/03/2010
///	Revision:	  ---
///=====================================================================================
///
int two_temperature_pulse(){
	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){std::cout << "program::two_temperature_pulse has been called" << std::endl;}
	
	// function prototypes

	// Setup LLG arrays
	sim::initialise();

	// Initialise spins to +z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		if(atoms::type_array[atom]==4){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=-1.0;
		}
		else{
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
		}		
	}
	vout::pov_file();
	// Set up loop variables
	sim::H_applied=0.0;
	
	const double Ce = 7.0E02; //electron specific heat 
	const double Cl = 3.0E06; //photon specific heat
	const double G = 17.0E17 ;//electron coupling constant
	
	double pump_time=20.0e-15; // Seconds
	double pump_power=2.4e22; // ?
	
	//for (int mat=0;mat<mp::num_materials;mat++){
	//for (int nmat=0;nmat<mp::num_materials;nmat++){
	//	std::cout << mp::material[mat].Jij_matrix[nmat] << std::endl;
	//}
	//}
	
	// Set initial system temperature
	sim::temperature=77.0;
	
	// Equilibrate system
	for(sim::time=0;sim::time<sim::equilibration_time;sim::time++){
		double actual_time = double(sim::time-sim::equilibration_time)*mp::dt_SI;
		LLG(1);
		
		// Calculate mag_m, mag
		if(sim::time%sim::partial_time==0){
			stats::mag_m();
			// Output to screen and file after each temperature
			if(vmpi::my_rank==0){
			  vmag << actual_time << "\t" << sim::temperature << "\t" << stats::total_mag_m_norm;
				vmag << "\t" << stats::total_mag_norm[0];
				vmag << "\t" << stats::total_mag_norm[1];
				vmag << "\t" << stats::total_mag_norm[2];
				//
				// loop over all materials and output sublattice moments
				for (int mat=0;mat<mp::num_materials;mat++){
					vmag << "\t" << stats::sublattice_mx_array[mat];
					vmag << "\t" << stats::sublattice_my_array[mat];
					vmag << "\t" << stats::sublattice_mz_array[mat];
					vmag << "\t" << stats::sublattice_magm_array[mat];
				}
				vmag << std::endl;
			}
		}
	}
	vout::pov_file();

	std::cout << "Equilibration complete" << std::endl;

	double Te=sim::temperature;
	double Tp=sim::temperature;

	//std::cout << "Timestep: "<< mp::dt_SI << std::endl;
	
	std::ofstream pp;
	pp.open("pp.dat");
	
	// Simulate system with single timestep resolution
	for(sim::time=0;sim::time<sim::loop_time;sim::time++){
		
		// calculate real time and temperature using gaussian cooling
		double actual_time = double(sim::time)*mp::dt_SI;
		//sim::temperature=RT+max_dT*exp(-((actual_time*actual_time)/(cooling_time*cooling_time)));
		double pump=pump_power*exp(-((actual_time-3.*pump_time)/(pump_time) )*((actual_time-3.*pump_time)/(pump_time) ));
		//double pump=pump_power*exp(-((actual_time*actual_time)/((pump_time-3*)*pump_time)));
		Te = (-G*(Te-Tp)+pump)*mp::dt_SI/(Ce*Te) + Te;
		Tp = ( G*(Te-Tp)     )*mp::dt_SI/Cl + Tp;
		pp << actual_time << "\t" << Te << "\t" << Tp << "\t" << pump << std::endl;
		
		sim::temperature=Te;
		// Calculate LLG
		LLG(1);
		
		// Calculate mag_m, mag
		if(sim::time%sim::partial_time==0){
			stats::mag_m();
			// Output to screen and file after each temperature
			if(vmpi::my_rank==0){
			  std::cout << actual_time << "\t";
			  std::cout << sim::temperature << "\t" << stats::total_mag_m_norm;
			  std::cout << "\t" << stats::total_mag_norm[0];
			  std::cout << "\t" << stats::total_mag_norm[1];
			  std::cout << "\t" << stats::total_mag_norm[2];
			  std::cout << std::endl;
			  vmag << actual_time << "\t" << sim::temperature << "\t" << stats::total_mag_m_norm;
				vmag << "\t" << stats::total_mag_norm[0];
				vmag << "\t" << stats::total_mag_norm[1];
				vmag << "\t" << stats::total_mag_norm[2];
				//
				// loop over all materials and output sublattice moments
				for (int mat=0;mat<mp::num_materials;mat++){
					vmag << "\t" << stats::sublattice_mx_array[mat];
					vmag << "\t" << stats::sublattice_my_array[mat];
					vmag << "\t" << stats::sublattice_mz_array[mat];
					vmag << "\t" << stats::sublattice_magm_array[mat];
				}
				vmag << std::endl;
			}
		}
		//if(sim::time%60==0){
		// vout::pov_file();
		// }
	}
	//vout::pov_file();
		vout::pov_file();
	return EXIT_SUCCESS;
}

/// @brief Function to calculate a static hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// accoring to the input flag - either randomly or ordered.For the ordered case the temperature sequence
/// increases from zero, for the random case the temperature decreases from the maximum temperature. After
/// initialisation the sytem is equilibrated for sim::equilibration timesteps.
///
/// @section notes Implementation Notes
/// Capable of hot>cold or cold>hot calculation. 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    28/01/2010
///
/// @param[in] init Determines whether the system is initialised randomly (0) or ordered (1)
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		28/01/2010
///	Revision:	  ---
///=====================================================================================
///
int static_hysteresis(){
	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){std::cout << "program::static_hysteresis has been called" << std::endl;}
	
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << " Running Static Hysteresis Loop Program" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	
	// Setup LLG arrays
	sim::initialise();

	// Initialise spins to near-ordered state
	for(int atom =0;atom<atoms::num_atoms;atom++){
			atoms::x_spin_array[atom]=0.0+0.05*mtrandom::grnd()-0.1;
			atoms::y_spin_array[atom]=0.0;
			atoms::z_spin_array[atom]=1.0;
	}
	
	// Disable temperature as this will prevent convergence
	sim::temperature = 0.0;
	sim::hamiltonian_simulation_flags[3] = 0;	// Thermal
	// Also set up high alpha to increase convergence rate
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].alpha=4.0;
	}
	// Set H vector
	sim::H_vec[0]=0.01;
	sim::H_vec[1]=0.0;
	sim::H_vec[2]=0.9999;
	
	double mod_H=1.0/sqrt(sim::H_vec[0]*sim::H_vec[0]+sim::H_vec[1]*sim::H_vec[1]+sim::H_vec[2]*sim::H_vec[2]);

	sim::H_vec[0]*=mod_H;
	sim::H_vec[1]*=mod_H;
	sim::H_vec[2]*=mod_H;
	
	// Estimate FM Corecivity
	double Hc = 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI;
	//std::cout << std::cout.flags() << std::endl;
	//std::cout.unsetf(std::ios_base::floatfield);
	//std::cout << std::cout.flags() << std::endl;
	
	std::cout << "\tEstimated Coercivity for FM is: " << Hc << " Tesla" << std::endl;
	std::cout << "\tmu_s: " << mp::material[0].mu_s_SI << std::endl;
	std::cout << "\tKu  : " << mp::material[0].Ku1_SI << std::endl;
	
	// Equilibrate system in strong positive field
	sim::H_applied=5.0;
	LLG(sim::equilibration_time);
	
	std::cout << "\t------------------------------------------------------------------------" << std::endl;
	std::cout << "\tEquilibration Complete" << std::endl;
	std::cout << "\t------------------------------------------------------------------------" << std::endl;
	
	// Output initial spin configuration
	vout::pov_file();
	
	// Setup min and max fields and increment (mT)
	int iHmax=round(sim::Hmax*1.0E3);
	int iHmin=round(sim::Hmin*1.0E3);
	int iHinc=round(sim::Hinc*1.0E3);
	// Perform Field Loop
	for(int parity=-1;parity<2;parity+=2){
		for(int H=iHmin;H<=iHmax;H+=iHinc){
			
			// Set applied field (Tesla)
			sim::H_applied=double(H)*double(parity)*1.0e-3;
			
			//sim::partial_time=1;
			// Simulate system
			for(sim::time=0;sim::time<sim::total_time;sim::time+=sim::partial_time){
				
				// Calculate LLG
				LLG_relax(sim::partial_time);
				

				//std::cout << atoms::x_spin_array[0] << "\t" << atoms::y_spin_array[0] << "\t" << atoms::z_spin_array[0] << std::endl;
			
				//std::cin.get();
				double torque=stats::max_torque();
				if((torque<1.0e-5) && (sim::time>100)){
					break;
				}

			}
			
			// Calculate mag_m, mag
			stats::mag_m();
			
			//std::cout << sim::time << "\t" << sim::H_applied << "\t" << atoms::x_spin_array[0] << "\t" << atoms::y_spin_array[0] << "\t" << atoms::z_spin_array[0] << "\t" << stats::max_torque() << std::endl;
			// Output to screen and file after each field
			if(vmpi::my_rank==0){
				std::cout << "\t" << sim::time << "\t" << sim::H_applied << "\t" << stats::total_mag_norm[0];
				std::cout << "\t" << stats::total_mag_norm[1] << "\t" << stats::total_mag_norm[2];
				std::cout << "\t" << stats::total_mag_m_norm << std::endl;

				vmag << sim::time << "\t" << sim::H_applied << "\t" << stats::total_mag_norm[0];
				vmag << "\t" << stats::total_mag_norm[1] << "\t" << stats::total_mag_norm[2];
				vmag << "\t" << stats::total_mag_m_norm << std::endl;
			}
			
		} // End of field loop
} // End of parity loop
	return EXIT_SUCCESS;
}

/// @brief Function to calculate the hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised 
/// ordered. After initialisation a whole hysteresis loop of the system and coercivity are calculated.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Weijia Fan, wf507@york.ac.uk
/// @version 1.0
/// @date    27/01/2010
///
/// @return EXIT_SUCCESS
/// 
/// @internal
///	Created:		27/01/2010
///	Revision:	  ---
///=====================================================================================
///


int hysteresis(){
	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){std::cout << "hysteresis has been called" << std::endl;}
	
	// Declare function prototype
	//int calculate_applied_fields(const int,const int);
	// set loop temperature
	// Setup LLG arrays
	sim::initialise();
	
	// Initialise spins to +z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
		//atoms::z_spin_array[atom]=1.0-2.0*double(vmpi::my_rank%2);
	}
	
	 // Setup Hmax and J=Hinc
	 int Hmax = 6000; // mT
	 int Hinc = 10; // mT	 
	 
	 // setup mag_new(H_new) and mag_old(H_old) for coercive field
	 double mag_new = 0.0;
	 double mag_old = 0.0;
	 double H_old = 0.0;
	 double H_new = 0.0;
	 double H_c_left= 0.0;
	 double H_c_right = 0.0;
	 double H_coercive = 0.0;
	 //std::cout << std::cout.flags() << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);
        //std::cout << std::cout.flags() << std::endl;
	 //std::cout << mp::material[0].mu_s_SI << "\t" << mp::material[0].Ku1_SI << "\t" << 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI << std::endl;
	 std::cout << "Estimated Coercivity:" << 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI << std::endl;
	 vout::pov_file();
	 // parity loop
	 for(int parity=-1;parity<2;parity+=2){
	  // Set up loop variables
	  for (int H = -Hmax; H<= Hmax;H+=Hinc){

	    sim::H_applied=double(parity)*double(H)*0.001;	// Tesla
	  
	  // time loop
	    for(sim::time=0;sim::time<sim::loop_time;sim::time+=sim::partial_time){
	      LLG(sim::partial_time);
	      stats::mag_m();
	      //if(vmpi::my_rank==0){
	      // 	std::cout << sim::time<< "\t" << stats::total_mag_m_norm;
	      //	std::cout << "\t" << stats::total_mag_norm[0];
	      //	std::cout << "\t" << stats::total_mag_norm[1];
	      //	std::cout << "\t" << stats::total_mag_norm[2];
	      //	std::cout << std::endl;
	      //}
	    } 
	    
		// output pov_file after each field point
	    //vout::pov_file();
		
	    std::cout << vmpi::my_rank;
	    std::cout << "\t" << stats::total_mag_m_norm;
	    std::cout << "\t" << stats::total_mag_norm[0];
	    std::cout << "\t" << stats::total_mag_norm[1];
	    std::cout << "\t" << stats::total_mag_norm[2];
	    std::cout << std::endl;
	    
	     if(vmpi::my_rank==0){
		mag_new = stats::total_mag_norm[2];
		H_new = sim::H_applied;
		if (mag_new*mag_old < 0 & mag_old > mag_new){
		  // calculate coercive field
		  H_c_left=(mag_old*H_new-mag_new*H_old)/(mag_old-mag_new);
		  std::cout << "\t" << "the left coercivity is" << "\t" << H_c_left << "\t" << "Tesla" << std::endl;	
		}
		if (mag_new*mag_old < 0 & mag_old < mag_new){
		  H_c_right=(mag_old*H_new-mag_new*H_old)/(mag_old-mag_new);
		  std::cout << "\t" << "the right coercivity is" << "\t" <<  H_c_right << "\t" << "Tesla" << std::endl;
		}		 		  
	      std::cout << sim::H_applied << "\t" << stats::total_mag_m_norm; // Tesla
	      std::cout << "\t" << stats::total_mag_norm[0];
	      std::cout << "\t" << stats::total_mag_norm[1];
	      std::cout << "\t" << stats::total_mag_norm[2];
	      std::cout << std::endl;
		// check current and before values
		//std::cout << "current mag_new is" << mag_new << "\t" << "before mag_old is" << mag_old<< std::endl;
		//std::cout << "current H_new is" << H_new << "\t" << "before H is" << H_old << std::endl;
	      vmag << sim::H_applied << "\t" << stats::total_mag_m_norm;
	      vmag << "\t" << stats::total_mag_norm[0];
	      vmag << "\t" << stats::total_mag_norm[1];
	      vmag << "\t" << stats::total_mag_norm[2];
	      vmag << std::endl;
	      

	      				
	      mag_old = mag_new;
	      H_old = H_new;
	    }
		 

	  //if ((H%100)==0){vout::pov_file();}		
	     // mag_new = stats::total_mag_norm[2];
	    
	     // if ((mag_old-0.8)/(mag_new-0.8)<0){
		//vout::pov_file();	      
	     // }
	     // if ((mag_old+0.8)/(mag_new+0.8)<0){
		//vout::pov_file();	      
	     // }
      	     // mag_old = mag_new;

	  }


	 }
	 if(vmpi::my_rank==0){
	  H_coercive = -H_c_left;   //0.5*(H_c_right-H_c_left);
	  std::cout << "coercive field of the system is" << "\t" << H_coercive << "\t" << "Tesla" << std::endl;
	  //vmag << "Hc+ is" << "\t" << H_c_right << "\tTesla" << "Hc- is" << "\t" << H_c_left << "\tTesla" << std::endl;
	  vmag << "The coercive field of the system is" << "\t" << H_coercive << "\t" << "Tesla" << std::endl;
	 }
	 
	return EXIT_SUCCESS;
  }
// End of namespace
}
