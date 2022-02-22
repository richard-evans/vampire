//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"

namespace program{
   void localised_temperature_pulse(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::localised_temperature_pulse has been called" << std::endl;}

   // Check localised temperature pulse is enabled
   if(!ltmp::is_enabled()){
      zlog << zTs() << "Fatal Error: Localised temperature pulse is not initialised. Exiting" << std::endl;
      std::cerr << "Fatal Error: Localised temperature pulse is not initialised. Exiting" << std::endl;
      err::vexit();
   }

   // Set equilibration temperature and field
   sim::temperature=sim::Teq;

   // Initialise electron and phonon temperature
   sim::TTTe=sim::temperature;
   sim::TTTp=sim::temperature;

   // Equilibrate system
   while(sim::time<sim::equilibration_time){

      // loop over partial_time to update temperature every time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate temperature
         ltmp::update_localised_temperature(-1.e-7);

         sim::integrate(1);

      }

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();
   }

   // record starting time after equilibration
   uint64_t start_time = sim::time;

   // Simulate temperature pulse
   while(sim::time < sim::total_time+start_time){

      // loop over partial_time to update temperature every time
      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         // Calculate time from pulse
         double time_from_start=mp::dt_SI*double(sim::time-start_time);

         // Calculate temperature
         ltmp::update_localised_temperature(time_from_start);

         // Integrate system
         sim::integrate(1);

      }

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

} // end of temperature pulse

} // end of namespace ltmp
