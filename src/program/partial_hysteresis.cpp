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

/// @brief Function to calculate the partial hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details sim:program=partial-hysteresis-loop simulates a partial hysteresis loop, starting at
///          sim:minimum-applied-field-strength to sim:maximum-applied-field-strength in steps of
///          sim:applied-field-strength-increment. Note that the sign of the increment is
///          significant, indicating the direction of the loop. Invalid combinations (which lead to
///          an infinite loop) are checked during the initialisation and will print out a warning.
///          As with the full hysteresis loop, the minimum resolution of the applied field
///          increment is 1 uT.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    02/12/20103
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		02/12/2013
///	Revision:	  ---
///=====================================================================================
///
void partial_hysteresis_loop(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::partial-hysteresis has been called" << std::endl;}

   // Equilibrate system in saturation field
   sim::H_applied=sim::Heq;
   sim::integrate(sim::equilibration_time);

   // Setup min and max fields and increment (uT)
   int64_t iHmax=vmath::iround64(double(sim::Hmax)*1.0E6);
   int64_t iHmin=vmath::iround64(double(sim::Hmin)*1.0E6);
   int64_t iHinc=vmath::iround64(double(sim::Hinc)*1.0E6);

   // Check for loop direction and adjust parameters
   // so that field loop works in positive sense
   double parity=1.0;
   if(iHinc < 0){
      iHmax=-iHmax;
      iHmin=-iHmin;
      iHinc=-iHinc;
      parity=-1.0;
   }

   // Perform Field Loop
   for(int H=iHmin;H<=iHmax;H+=iHinc){

      // Set applied field (Tesla)
      sim::H_applied=double(H)*parity*1.0e-6;

      // Reset start time
      uint64_t start_time=sim::time;

      // Reset mean magnetisation counters
      stats::reset();

      // Integrate system
      while(sim::time<sim::loop_time+start_time){

         // Integrate system
         sim::integrate(sim::partial_time);

         // Calculate mag_m, mag
         stats::update();

      }

      // Output to screen and file after each field
      vout::data();

   } // End of field loop

   return;
}

}//end of namespace program
