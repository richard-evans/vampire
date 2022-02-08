//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cstdlib>

// Vampire Header files
#include "vmath.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"

namespace program{

/// @brief Function to calculate a static hysteresis loop
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk, Andrea Meo, am1808@york.ac.uk (Revision)
/// @version 1.0
/// @date    10/03/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		28/01/2010
///	Revision:	06/07/2017
///=====================================================================================
///
int static_hysteresis(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::static_hysteresis has been called" << std::endl;}

	// Disable temperature as this will prevent convergence
	sim::temperature = 0.0;
	sim::hamiltonian_simulation_flags[3] = 0;	// Thermal

	// Setup min and max fields and increment (uT)
	int64_t iHmax=vmath::iround64(double(sim::Hmax)*1.0E6);
	int64_t miHmax=-iHmax;
	int64_t parity_old;
	int64_t iH_old;
	int64_t start_time;

	// Equilibrate system in saturation field, i.e. the largest between equilibration and maximum field set by the user
   if(sim::Heq >= sim::Hmax){
	   sim::H_applied=sim::Heq;
   }
   else{
   	sim::H_applied=sim::Hmax;
   }

	// Initialise sim::integrate only if it not a checkpoint
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else sim::integrate(sim::equilibration_time);

   // Hinc must be positive
	int64_t iHinc=vmath::iround64(double(fabs(sim::Hinc))*1.0E6);

   int64_t Hfield;
   int64_t iparity=sim::parity;
	parity_old=iparity;

   // Save value of iH from previous simulation
	if(sim::load_checkpoint_continue_flag) iH_old=int(sim::iH);

	// Perform Field Loop -parity
	while(iparity<2){
      // If checkpoint is loaded with continue flag, then set up correctly max,min field values
		if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag)
      {
         //necessary to upload value of iH_old when loading the checkpoint !!!
		   iH_old=static_cast<int64_t>(sim::iH); //  int(sim::iH);
         //Setup min and max fields and increment (uT)
			if(parity_old<0){
				if(iparity<0) miHmax=iH_old;
				else if(iparity>0 && iH_old<=0) miHmax=iH_old; //miHmax=(iHmax-iHinc);
				else if(iparity>0 && iH_old>0) miHmax=-(iHmax);
			}
			else if(parity_old>0) miHmax=iH_old;
			Hfield=miHmax;
		}
		else	Hfield=miHmax;

		// Perform Field Loop -field
		while(Hfield<=iHmax){

			// Set applied field (Tesla)
			sim::H_applied=double(Hfield)*double(iparity)*1.0e-6;

			// Reset start time
			start_time=sim::time;

			// Reset mean magnetisation counters
			stats::mag_m_reset();

			// Integrate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);

				double torque=stats::max_torque(); // needs correcting for new integrators
				if((torque<1.0e-6) && (sim::time-start_time>100)){
					break;
				}

				// Calculate mag_m, mag
				stats::mag_m();

         } // End integration loop

			// Increment of iH
			Hfield+=iHinc;
			sim::iH=Hfield; //sim::iH+=iHinc;


			// Output to screen and file after each field
			vout::data();

		} // End of field loop

		// Increment of parity
      iparity+=2;
      sim::parity=iparity;

	} // End of parity loop

	return EXIT_SUCCESS;
} // end of static-hysteresis-loop

}//end of namespace program
