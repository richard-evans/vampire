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
/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@nanohpc.com
/// @version 1.0
/// @date    10/06/2011
/// @internal
///   Created:    10/06/2011
///   Revision:   Andrea Meo, Feb 2022
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "hamr.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
/// @details Performs a time series with moving localised head field and temperature pulse
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    10/06/2011
///
/// @internal
///   Created:    10/06/2011
///   Revision:   Andrea Meo, 2022
///=====================================================================================
///
   void hamr(){

      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "program::hamr has been called" << std::endl;}

         // Set partial_time (sim:time-steps-increment) to 1 to ensure good simulation
         sim::partial_time = 1;

         // Set equilibration temperature and field
         sim::temperature=sim::Teq;

         // Disable laser
         hamr::head_laser_on=false;

         // Equilibrate system
         while(sim::time<sim::equilibration_time){

            sim::integrate(sim::partial_time);

            // Calculate magnetisation statistics
            stats::update();

            // Output data
            vout::data();
         }

         // now enable laser
         hamr::head_laser_on=true;
         std::cout << " Equilibration terminated, enabling laser ..." << std::endl;
         zlog << zTs() << "Equilibration terminated, enabling laser ..." << std::endl;

         // Set system temperature as minimum temperature
         sim::temperature=sim::Tmin;

         // Reset mean magnetisation counters (unless continuing from a checkpoint saved after the reset)
         if(!(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::time > sim::equilibration_time)) stats::reset();

         // Perform harm continuous simulation
         hamr::hamr_continuous();

         // //-------------------------------------------------------------------------//
         // // force outputting at end of simulation after laser has been switched off
         // //-------------------------------------------------------------------------//
         // // Disable laser
         // hamr::head_laser_on=false;
         // // Switch off external field
         // sim::H_applied = 0.0;
         // // Set system temperature as minimum temperature
         // sim::temperature=sim::Tmin;
         // std::cout << " Disabling laser and external field and integrating system for 1 time-step" << std::endl;
         // zlog << zTs() << " Disabling laser and external field and integrating system for 1 time-step" << std::endl;
         // // Integrate
         // sim::integrate(sim::partial_time);
         // // Calculate magnetisation statistics
         // stats::update();
         // // Output data
         // vout::data();
         // std::cout << " Outputting system at the end of HAMR simulations\n" << std::endl;
         // zlog << zTs() << "Outputting system at the end of HAMR simulations" << std::endl;

   } // end of hamr()

}//end of namespace program
