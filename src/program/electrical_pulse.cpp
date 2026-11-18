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
#include <cmath>

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
// Square pulse: linear rise, constant strength for electrical-pulse-time, linear
// fall. Time t is measured from the start of the pulse.
//------------------------------------------------------------------------------
static double electrical_pulse_square(const double t){

   const double rise = pgi::electrical_pulse_rise_time;
   const double on   = pgi::electrical_pulse_time;
   const double fall = pgi::electrical_pulse_fall_time;

   double fraction = 0.0;

   if( t < 0.0 )                  fraction = 0.0;                              // before pulse
   else if( t < rise )            fraction = t / rise;                         // rise time
   else if( t < rise + on )       fraction = 1.0;                              // continuous current
   else if( t < rise + on + fall) fraction = 1.0 - (t - (rise + on)) / fall;   // fall time
   else                           fraction = 0.0;                              // after pulse

   return pgi::electrical_pulse_strength * fraction;

}

namespace program{
namespace internal{

//------------------------------------------------------------------------------
// Gaussian pulse with full width at half maximum electrical-pulse-time, peaking at
// t = 3 x electrical-pulse-time. This matches the laser Gaussian in
// two_temperature_function, so equal pulse times and zero delay overlap the
// laser and electrical pulses. Time t is measured from the end of the delay.
//------------------------------------------------------------------------------
double electrical_pulse_gaussian(const double t){

   const double tau      = electrical_pulse_time;
   const double four_ln2 = 2.77258872224; // 4 ln 2
   const double reduced  = (t - 3.0 * tau) / tau;

   return electrical_pulse_strength * std::exp(-four_ln2 * reduced * reduced);

}

//------------------------------------------------------------------------------
// Fractional electric field strength (scales STT/SOT and voltage) at a time from
// the start of the pulse program: zero during the delay, then the selected pulse
// shape, repeated every electrical-pulse-period if one is set.
//------------------------------------------------------------------------------
double electrical_pulse_fraction(const double time_from_start){

   // shift time origin by user-specified delay (zero by default)
   double t = time_from_start - electrical_pulse_delay;

   // before pulse delay: current = 0
   if( t < 0.0 ) return 0.0;

   // repeated pulses
   if( electrical_pulse_period > 0.0 ) t = std::fmod(t, electrical_pulse_period);

   switch(electrical_pulse_shape){
      case electrical_pulse_shape_t::square:
         return electrical_pulse_square(t);
      case electrical_pulse_shape_t::gaussian:
         return electrical_pulse_gaussian(t);
   }

   return 0.0;

}

} // end of namespace internal

//------------------------------------------------------------------------------
// Electrical pulse program
//------------------------------------------------------------------------------
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
         pg::fractional_electric_field_strength = pgi::electrical_pulse_fraction(time_from_start);

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
