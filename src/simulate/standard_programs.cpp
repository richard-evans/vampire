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
#include "errors.hpp"
#include "grains.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

int initialise_system();
	
/// @namespace program
/// @brief A Namespace containing functions for standard programs.
/// 
/// @internal
///=====================================================================================
///
namespace program{








/*/// @brief Function to calculate a time series with gaussian cooling profile
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
	if(err::check==true){std::cout << "program::hamr_run has been called" << std::endl;}
	
	// function prototypes

        // perform 10 sequential runs
        for(int run=0;run<1;run++){

	// Initialise spins to +z
	for(int atom =0;atom<atoms::num_atoms;atom++){
		atoms::x_spin_array[atom]=0.0;
		atoms::y_spin_array[atom]=0.0;
		atoms::z_spin_array[atom]=1.0;
	}
	
	// Set up loop variables
	sim::H_applied=-0.8;
	
	double cooling_time=0.5e-9; // Seconds
	double max_dT = sim::delta_temperature; 
	double RT = 300.0;
	// Set initial system temperature
	sim::temperature=300.0;

        // Initialise random number generator
        mtrandom::grnd.seed(vmpi::my_rank);

	// Equilibrate system
	//sim::LLG(sim::equilibration_time);

	// Simulate system with single timestep resolution
	for(sim::time=0;sim::time<sim::loop_time;sim::time++){
		
		// calculate real time and temperature using gaussian cooling
		double actual_time = (double(run)*double(sim::loop_time)+double(sim::time))*mp::dt_SI;
		double rtime = double(sim::time)*mp::dt_SI;
		sim::temperature=RT+max_dT*exp(-((rtime*rtime)/(cooling_time*cooling_time)));
		
		// Calculate LLG
		//sim::LLG(1);
		
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
        }
	//vout::pov_file();
	// output final grain magnetisations
	std::ofstream vgrain;
	vgrain.open("vgrain");
	grains::output_mag(vgrain);
	vgrain.close();
		
	return EXIT_SUCCESS;
}*/
/*
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
	if(err::check==true){std::cout << "program::two_temperature_pulse has been called" << std::endl;}
	
	// function prototypes

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
	//vout::pov_file();
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
		//sim::LLG(1);
		
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
	//vout::pov_file();

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
		//sim::LLG(1);
		
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
		//vout::pov_file();
	return EXIT_SUCCESS;
}*/

}//end of namespace program



