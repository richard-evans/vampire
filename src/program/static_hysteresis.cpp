///
/// @file
/// @brief Contains the Static Hysteresis program
///
/// @details Performs a field loop to determine field dependent magnetic behaviour
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    10/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	10/03/2011
///=====================================================================================
///

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// Standard Libraries
#include <iostream>
#include <cmath>

namespace program{
	
/// @brief Function to calculate a static hysteresis loop
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    10/03/2011
///
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
	
	// Disable temperature as this will prevent convergence
	sim::temperature = 0.0;
	sim::hamiltonian_simulation_flags[3] = 0;	// Thermal

	// Also set up high alpha to increase convergence rate
	for(int mat=0;mat<mp::num_materials;mat++){
		mp::material[mat].alpha=4.0;
	}

	// Ensure H vector is unit length
	double mod_H=1.0/sqrt(sim::H_vec[0]*sim::H_vec[0]+sim::H_vec[1]*sim::H_vec[1]+sim::H_vec[2]*sim::H_vec[2]);
	sim::H_vec[0]*=mod_H;
	sim::H_vec[1]*=mod_H;
	sim::H_vec[2]*=mod_H;
	
	// Estimate Material Coercivities
	double Hc = 2.0*mp::material[0].Ku1_SI/mp::material[0].mu_s_SI;
	
	std::cout << "\tEstimated Coercivity for FM is: " << Hc << " Tesla" << std::endl;
	std::cout << "\tmu_s: " << mp::material[0].mu_s_SI << std::endl;
	std::cout << "\tKu  : " << mp::material[0].Ku1_SI << std::endl;
	
	// Equilibrate system in strong positive field
	sim::H_applied=sim::Hmax;
	sim::LLG(sim::equilibration_time);
		
	// Setup min and max fields and increment (mT)
	int iHmax=iround(double(sim::Hmax)*1.0E3);
	int iHmin=iround(double(sim::Hmin)*1.0E3);
	int iHinc=iround(double(sim::Hinc)*1.0E3);

	// Perform Field Loop
	for(int parity=-1;parity<2;parity+=2){
		for(int H=iHmin;H<=iHmax;H+=iHinc){
			
			// Set applied field (Tesla)
			sim::H_applied=double(H)*double(parity)*1.0e-3;
			
			// Simulate system
			while(sim::time<sim::loop_time){
			
				// Integrate system
				sim::integrate(sim::partial_time);
				
				double torque=stats::max_torque();
				if((torque<1.0e-5) && (sim::time>100)){
					break;
				}

			}
			
			// Calculate mag_m, mag
			stats::mag_m();
			
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

