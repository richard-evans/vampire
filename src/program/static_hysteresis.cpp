#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include <iostream>

namespace program{
	
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
	if(err::check==true){std::cout << "program::static_hysteresis has been called" << std::endl;}
	
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
	sim::LLG(sim::equilibration_time);
	
	std::cout << "\t------------------------------------------------------------------------" << std::endl;
	std::cout << "\tEquilibration Complete" << std::endl;
	std::cout << "\t------------------------------------------------------------------------" << std::endl;
	
	// Output initial spin configuration
	vout::pov_file();
	
	// Setup min and max fields and increment (mT)
	int iHmax=round(double(sim::Hmax)*1.0E3);
	int iHmin=round(double(sim::Hmin)*1.0E3);
	int iHinc=round(double(sim::Hinc)*1.0E3);
	// Perform Field Loop
	for(int parity=-1;parity<2;parity+=2){
		for(int H=iHmin;H<=iHmax;H+=iHinc){
			
			// Set applied field (Tesla)
			sim::H_applied=double(H)*double(parity)*1.0e-3;
			
			//sim::partial_time=1;
			// Simulate system
			for(sim::time=0;sim::time<sim::total_time;sim::time+=sim::partial_time){
				
				// Calculate LLG
				sim::LLG_relax(sim::partial_time);
				

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


}//end of namespace program

