//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2022. All rights reserved.
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

// Vampire headers
//
//#include "ltmp.hpp"
//#include "material.hpp"


// namespace abbreviation for brevity
namespace pg = program;
namespace pgi = program::internal;

//------------------------------------------------------------------------------
// Function to update STT/SOT and spin transport parameters
//------------------------------------------------------------------------------
void update_electric_field_strength(const double time_from_start){

   // implement rise time
   if( time_from_start < pgi::electrical_pulse_rise_time ){
      pg::fractional_electric_field_strength = time_from_start / pgi::electrical_pulse_rise_time;
   }
   // implement continuous current
   else if( time_from_start < pgi::electrical_pulse_rise_time + pgi::electrical_pulse_time ){
      pg::fractional_electric_field_strength = 1.0;
   }
   // implement fall time
   else if( time_from_start < pgi::electrical_pulse_rise_time + pgi::electrical_pulse_time + pgi::electrical_pulse_fall_time) {
      const double fractional_fall_time = time_from_start - (pgi::electrical_pulse_rise_time + pgi::electrical_pulse_time);
      pg::fractional_electric_field_strength = 1.0 - fractional_fall_time / pgi::electrical_pulse_fall_time;
   }
   // after pulse current = 0
   else{
      pg::fractional_electric_field_strength = 0.0;
   }

   return;

}


namespace program{


void electrical_pulse(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::electrical_pulse has been called" << std::endl;}

   // Set equilibration temperature and zero current
   const double temp = sim::temperature; // current simulation temperature
   sim::temperature = sim::Teq;
   pg::fractional_electric_field_strength = 0.0;

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

   // Simulate electrical pulse
   while(sim::time < sim::total_time+start_time){

      // loop over partial_time to update temperature every time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate time from pulse
         double time_from_start = mp::dt_SI * double(sim::time-start_time);

         // Calculate electrical field strength
         update_electric_field_strength(time_from_start);

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
