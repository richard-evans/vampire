//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
//#include "atoms.hpp"
//#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
//#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
//#include "vmath.hpp"
//#include "vmpi.hpp"

namespace program{

//----------------------------------------------------------------------------------------------
// Function to calculate field cooling with different maximum and minimum
// temperatures for each material
//----------------------------------------------------------------------------------------------
void local_field_cool(){

	// Perform several runs if desired
	for(int run=0;run<sim::runs; run++){

		// Set equilibration temperature and field
		sim::temperature=sim::Teq;

		// Equilibrate system
		while(sim::time<sim::equilibration_time){

			sim::integrate(sim::partial_time);

			// Calculate magnetisation statistics
			stats::update();

			// Output data
			vout::data();
		}

		uint64_t start_time = sim::time;

		// Perform Field Cooling
		while(sim::time<sim::total_time+start_time){

         // loop over partial time
         for(uint64_t tt=0; tt < sim::partial_time; tt++){

            // Calculate normalised temperature using desired cooling function
            double time_from_start=mp::dt_SI*double(sim::time-start_time);
				double normalised_temperature = 0.0;

            switch(sim::cooling_function_flag){
               case 0:{ // exponential
                  normalised_temperature = 0.0 + exp(-time_from_start/sim::cooling_time);
               }
               break;
               case 1:{ // gaussian
                  normalised_temperature = 0.0 + exp(-(time_from_start)*(time_from_start)/((sim::cooling_time)*(sim::cooling_time)));
               }
               break;
               case 2:{ // double-gaussian
                  double centre_time = 3.0*sim::cooling_time;
                  double time_from_centre=mp::dt_SI*double(sim::time-start_time)-centre_time;
                  normalised_temperature = 0.0 + exp(-(time_from_centre)*(time_from_centre)/((sim::cooling_time)*(sim::cooling_time)));
               }
               break;
               case 3:{ // linear
                  double Tlinear = 1.0-(time_from_start/sim::cooling_time);
                  if(Tlinear >= 0.0) normalised_temperature = Tlinear;
                  else normalised_temperature = 0.0;
               }
               break;
            }

				// set localised temperatures for each material
				for(int m = 0; m < mp::num_materials; ++m){
					mp::material[m].temperature = mp::material[m].minimum_temperature + normalised_temperature * (mp::material[m].maximum_temperature - mp::material[m].minimum_temperature);
				}

            // Integrate system
            sim::integrate(1);
         }

			// Calculate magnetisation statistics
			stats::update();

			// Output data
			vout::data();

		}

	} // end of run loop

} // end of field_cool()

}//end of namespace program
