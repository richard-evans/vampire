//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2015. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Time Series program
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///   Created:    30/03/2011
///   Revision:   --
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

namespace program{

/// @brief Function to calculate magnetisation over a time series
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
///
/// @internal
///   Created:    30/03/2011
///   Revision:   --
///=====================================================================================
///
void fmr(){

   // check calling of routine if error checking is activated
   if(err::check==true) std::cout << "program::time_series has been called" << std::endl;

   double temp=sim::temperature;

   // disable fmr fields for equilibration
   sim::enable_fmr = false;

   // Set equilibration temperature only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
      sim::temperature=sim::Teq;
   }

   // Equilibrate system
   while(sim::time<sim::equilibration_time){

      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();
   }

   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set simulation temperature
      sim::temperature=temp;

      // Reset mean magnetisation counters
      stats::reset();

   }

   // enable fmr fields
   sim::enable_fmr = true;

   // Initialize direction along z if not already set
   if(sim::fmr_field_unit_vector.size() != 3){
      sim::fmr_field_unit_vector.resize(3,0.0);
      sim::fmr_field_unit_vector[2] = 1.0;
   }

   // Perform Time Series
   while(sim::time<sim::equilibration_time+sim::total_time){

      // Integrate system
      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

}

}//end of namespace program
