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
/// @brief Contains the Temperature Pulse program
///
/// @details Performs a time series with variable (pulsed) temperature
///
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 2.0
/// @date    08/11/2012
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

/// function forward declarations
double two_temperature_function(double ftime);
double double_pump_two_temperature_function(double ftime);
double square_temperature_function(double ftime);
double double_pump_square_temperature_function(double ftime);

/// Selects appropriate temperature pulse function according to user input
double temperature_pulse_function(double function_time){

	switch(sim::pump_function){

		case(square): 
			return square_temperature_function(function_time);
			break;
		
		case(two_temperature):
			return two_temperature_function(function_time);
			break;
	
		case(double_pump_two_temperature):
			return double_pump_two_temperature_function(function_time);
			break;

		case(double_pump_square):
			return double_pump_square_temperature_function(function_time);
			break;

		default:
			return 0.0;
	}
}

///-----------------------------------------------------------------------------------------
///     Calculates temperature for single Gaussian pulse using the two-temperature model
///
///     Uses the free electron approximation in calculating the electron temperature
///     where C_e = gamma_e T_e (ref Unai Atxitia PhD thesis pp 53).
///     gamma_e ~ 3e3 J/m^-3/K^-2
///
///-----------------------------------------------------------------------------------------
double two_temperature_function(double ftime){

		const double reduced_time =  (ftime-3.*sim::pump_time)/(sim::pump_time);
		const double pump=sim::pump_power*exp(-reduced_time*reduced_time);

		const double Te = sim::TTTe;
		const double Tp = sim::TTTp;
		const double G  = sim::TTG;
		const double Ce = sim::TTCe;
		const double Cl = sim::TTCl;
		const double dt = mp::dt_SI;

		sim::TTTe = (-G*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
		sim::TTTp = ( G*(Te-Tp)     )*dt/Cl + Tp - (Tp-sim::Teq)*sim::HeatSinkCouplingConstant*dt;

      // Optionally set material specific temperatures
      if(sim::local_temperature==true){
         for(int mat=0;mat<mp::material.size();mat++){
            if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
            else mp::material[mat].temperature=sim::TTTe;
         }
      }

		return sim::TTTe;

}

/// Calculates temperature for double Gaussian pulses using the two-temperature model
double double_pump_two_temperature_function(double ftime){
	
		const double reduced_time =  (ftime-3.*sim::pump_time)/(sim::pump_time);
		const double pump1=sim::pump_power*exp(-reduced_time*reduced_time);

		const double reduced_time2 =  (ftime-sim::double_pump_delay-3.*sim::double_pump_time)/(sim::double_pump_time);
		const double pump2=sim::double_pump_power*exp(-reduced_time2*reduced_time2);

		const double Te = sim::TTTe;
		const double Tp = sim::TTTp;
		const double G  = sim::TTG;
		const double Ce = sim::TTCe;
		const double Cl = sim::TTCl;
		const double dt = mp::dt_SI;

		sim::TTTe = (-G*(Te-Tp)+pump1+pump2)*dt/(Ce*Te) + Te;
		sim::TTTp = ( G*(Te-Tp)           )*dt/Cl + Tp - (Tp-sim::Teq)*sim::HeatSinkCouplingConstant*dt;

      // Optionally set material specific temperatures
      if(sim::local_temperature==true){
         for(int mat=0;mat<mp::material.size();mat++){
            if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
            else mp::material[mat].temperature=sim::TTTe;
         }
      }

		return sim::TTTe;
}

/// Calculates temperature for single square pulse
double square_temperature_function(double ftime){

	//
	///
	///
	///            Tmax   |---------------|
	///                   |               |
	///                   |               |
	/// ------------------|               |--------------- Teq
	///                    <- pump time ->
	// 
	
	if(ftime <= sim::pump_time) return sim::Tmax;
	else return sim::Teq;
	
}

/// Calculates temperature for double square pulse
double double_pump_square_temperature_function(double ftime){

	///
	///          0               t1                t2              t3
	///
	///   Tmax   |---------------|
	///          |               |                 |---------------| double_pump_Tmax
	///          |               |                 |               |
	/// ---------|  -  -  -  -   |-----------------|   -   -   -   |------------ Teq
	///
	///           <- pump time -> <- dpump_delay -> <- pump time2 ->  
	///
	
	const double t1 = sim::pump_time;
	const double t2 = t1+sim::double_pump_delay;
	const double t3 = t2+sim::double_pump_time;
	
	if(ftime <= t1) return sim::Tmax;
	else if( (ftime > t2) && (ftime <=t3) ) return sim::double_pump_Tmax;
	else return sim::Teq;
	
}

namespace program{
/// @brief Function to calculate a time series with two temperature model heating/cooling profile
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    06/04/2011
/// 
/// @internal
///	Created:		10/03/2010
///	Revision:	06/04/2011
///=====================================================================================
///
void temperature_pulse(){
	
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::two_temperature_pulse has been called" << std::endl;}
		
	// Set equilibration temperature and field
	sim::temperature=sim::Teq;

	// Initialise electron and phonon temperature
	sim::TTTe=sim::temperature;
	sim::TTTp=sim::temperature;

   // If local temperature is set then also initalise local temperatures
   if(sim::local_temperature==true){
      for(int mat=0;mat<mp::material.size();mat++){
         mp::material[mat].temperature=sim::TTTe;
      }
   }

   // Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();
	}

	//loop sim::runs times
	for(int r=0; r<sim::runs;r++){
		
	// record starting time after equiibration/last pulse
	int start_time=sim::time;

	// Simulate temperature pulse
	while(sim::time<sim::total_time+start_time){

		// loop over partial_time to update temperature every time 
		for(int tt=0; tt < sim::partial_time; tt++){

			// Calculate time from pulse
			double time_from_start=mp::dt_SI*double(sim::time-start_time);
			
			// Calculate temperature
			sim::temperature=temperature_pulse_function(time_from_start);
			
			// Integrate system
			sim::integrate(1);
			
		}
		
		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();

	}
	
	}

} // end of temperature pulse

} // end of namespace program

