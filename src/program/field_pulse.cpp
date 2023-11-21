//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2023. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"

// program module headers
#include "internal.hpp"

// namespace abbreviation for brevity
namespace pg = program;
namespace pgi = program::internal;

//------------------------------------------------------------------------------
// Function to calculate a time-dependent field pulse
//------------------------------------------------------------------------------
namespace program{

void field_pulse(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::field_pulse has been called" << std::endl;}

   // Set equilibration temperature and zero field
   const double temp = sim::temperature; // current simulation temperature
   sim::temperature = sim::Teq;

   // Save input applied field strength
   const double max_field = sim::H_applied;

   // Equilibrate system
   while( sim::time < sim::equilibration_time){

      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

   // record starting time after equilibration
   uint64_t start_time = sim::time;

   // Set constant temperature
   sim::temperature = temp;

   // set centre time of field pulse
   const double centre_time = 3.0 * pgi::field_pulse_time;

   // set pulse_time^2
   const double pulse_time_sq = pgi::field_pulse_time * pgi::field_pulse_time;

   // Simulate field pulse
   while(sim::time < sim::total_time+start_time){

      // loop over partial_time to update temperature every time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate time from pulse
         double time_from_start = mp::dt_SI * double(sim::time-start_time);

         // Calculate applied field strength
         double time_from_centre = time_from_start - centre_time;
         sim::H_applied = max_field * exp(-(time_from_centre)*(time_from_centre)/pulse_time_sq);

         // Integrate system
         sim::integrate(1);

      }

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

   return;

} // end of electrical pulse

} // end of namespace program
