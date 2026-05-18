//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Implements the laser-electrical-pulse program, which combines a laser pulse
// (via the two-temperature model) with an independent Gaussian electrical pulse
// delayed by a user-specified time. The two pulses share no coupling — they are
// parameterised independently and may be shifted relative to each other.
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// program module headers
#include "internal.hpp"

namespace pg  = program;
namespace pgi = program::internal;

// temperature_pulse_function is defined in temperature_pulse.cpp and dispatches
// to the appropriate TTM or square profile selected by sim:laser-pulse-temporal-profile
extern double temperature_pulse_function(double function_time);

//------------------------------------------------------------------------------
// Gaussian envelope for the electrical pulse.
//
// The pulse is centred at (delay + 3*tau) from simulation start, matching the
// convention used by the laser Gaussian in two_temperature_function where the
// peak sits at t = 3*pump_time.  This means setting delay = 0 and
// electrical-pulse-time = laser-pulse-time produces perfectly overlapping pulses.
//------------------------------------------------------------------------------
static void update_laser_electrical_gaussian(const double t_elec){

   const double tau      = pgi::electrical_pulse_time;
   const double four_ln2 = 2.77258872224; // 4 ln 2
   const double reduced  = (t_elec - 3.0 * tau) / tau;
   pg::fractional_electric_field_strength = std::exp(-four_ln2 * reduced * reduced);

}

namespace program{

void laser_electrical_pulse(){

   if(err::check) std::cout << "program::laser_electrical_pulse has been called" << std::endl;

   // Initialise TTM temperatures at equilibrium and zero the electrical field
   sim::temperature = sim::Teq;
   sim::TTTe        = sim::Teq;
   sim::TTTp        = sim::Teq;
   pg::fractional_electric_field_strength = 0.0;

   // Initialise per-material temperatures when local temperature mode is active
   if(sim::local_temperature){
      for(unsigned int mat = 0; mat < mp::material.size(); mat++){
         if(mp::material[mat].couple_to_phonon_temperature) mp::material[mat].temperature = sim::TTTp;
         else                                               mp::material[mat].temperature = sim::TTTe;
      }
   }

   // Equilibrate system at Teq with no electrical field
   while(sim::time < sim::equilibration_time){
      sim::integrate(sim::partial_time);
      stats::update();
      vout::data();
   }

   const uint64_t start_time = sim::time;

   // Simulate combined laser + electrical pulse
   while(sim::time < sim::total_time + start_time){

      for(uint64_t tt = 0; tt < sim::partial_time; tt++){

         const double t = mp::dt_SI * double(sim::time - start_time);

         // Update spin temperature via TTM laser pulse
         sim::temperature = temperature_pulse_function(t);

         // Update Gaussian electrical pulse, offset by the user-specified delay
         update_laser_electrical_gaussian(t - pgi::electrical_pulse_delay);

         sim::integrate(1);
      }

      stats::update();
      vout::data();
   }

} // end of laser_electrical_pulse

} // end of namespace program
