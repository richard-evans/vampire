//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the standard benchmark program
///
/// @details Simulates a system for a number of timesteps
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///   Created:    05/02/2011
///   Revision:   09/03/2011
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace program{

int bmark(){
   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::bmark has been called" << std::endl;}

   // Calculate magnetisation statistics
   stats::update();

   // Output data
   vout::data();

   // Simulate system
   while(sim::time<sim::total_time){
      sim::integrate(sim::partial_time);

      // Calculate mag_m, mag after sim::partial_time steps
      stats::update();

      vout::data();

   } // end of time loop

   return EXIT_SUCCESS;
}

}//end of namespace program
