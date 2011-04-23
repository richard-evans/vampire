///
/// @file
/// @brief Contains the Hysteresis program
///
/// @details Performs a field loop to calculate field dependent magnetisation
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    21/03/2011
/// @internal
///	Created:		05/02/2011
///	Revision:	21/03/2011
///=====================================================================================
///

// Standard Libraries
#include <cstdlib>

// Vampire Header files
#include "vmath.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"


namespace program{

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
	if(err::check==true){std::cout << "program::hysteresis has been called" << std::endl;}
	
	// Equilibrate system in saturation field
	sim::H_applied=sim::Hmax;
	sim::integrate(sim::equilibration_time);
		
	// Setup min and max fields and increment (mT)
	int iHmax=vmath::iround(double(sim::Hmax)*1.0E3);
	int iHmin=vmath::iround(double(sim::Hmin)*1.0E3);
	int iHinc=vmath::iround(double(sim::Hinc)*1.0E3);

	// Perform Field Loop
	for(int parity=-1;parity<2;parity+=2){
		for(int H=iHmin;H<=iHmax;H+=iHinc){
			
			// Set applied field (Tesla)
			sim::H_applied=double(H)*double(parity)*1.0e-3;
			
			// Reset start time
			int start_time=sim::time;
			
			// Integrate system
			while(sim::time<sim::loop_time+start_time) sim::integrate(sim::partial_time);
			
			// Calculate mag_m, mag
			stats::mag_m();
			
			// Output to screen and file after each field
			vout::data();
			
		} // End of field loop
	} // End of parity loop

	return EXIT_SUCCESS;
  }

}//end of namespace program


