//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <math.h>   // for round() function

// Vampire headers
#include "errors.hpp"
#include "hamr.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"


// hamr headers
#include "internal.hpp"

namespace hamr{

	/* ------------------------------------------------------------------  /
	/  Continuous HAMR process                                             /
	/  T(x,y,t) = (Tmax-Tmin)* exp( -(x-v*t)*(x-v*t)/(2*sigmax*sigmax) )   /
	/                        * exp( -(y-Yc)*(y-Yc)  /(2*sigmay*sigmay) )   /
	/                                                                      /
	/  v=Head speed; Yc=Head y-coordinate                                  /
	/  sigmax,sigmay = FWHM in x,y-directions                              /
	/ ------------------------------------------------------------------- */
   void hamr_continuous(){

		// std::cout << "Performing HAMR continuous simulation" << std::endl;
		// zlog << zTs() << "Performing HAMR continuous simulation" << std::endl;
		
		// Determine times in terms of integration steps
		const double BL = hamr::internal::bit_size;
		const double TW = hamr::internal::track_size;
		const double speed = hamr::internal::head_speed;
		const double NPS = hamr::internal::NPS;
		const double padding = hamr::internal::track_padding;
		const double head_position_initial = 0.0; // -hamr::internal::system_dimensions_x*0.5 - NPS;
		const int n_bits = hamr::internal::num_bits;
		const int n_bits_per_tack = hamr::internal::bits_per_tack;
		const uint64_t bit_time   = int(round((BL/speed)/mp::dt_SI));
		const uint64_t extra_time = 0; //int(((BL*0.5)/speed)/mp::dt_SI);  // time to sweeping across half of bit size necessary to cover initial distance of head from track
		const uint64_t NPS_time   = int(round((NPS/speed)/mp::dt_SI)); // time to cover extra sweeping due to NPS
		const uint64_t track_time = n_bits_per_tack*bit_time + extra_time*2; 
		const uint64_t total_time =  n_bits * bit_time + extra_time*2;
		const uint64_t ramp_time = int(round(hamr::internal::H_ramp_time/mp::dt_SI));
		// Initial head position
		hamr::internal::head_position_x = head_position_initial;	
		hamr::internal::head_position_y = 0.5*hamr::internal::system_dimensions_y;
		// Initialise field magnitude to min value at beginning of simulation
		sim::H_applied = hamr::internal::Hmin;
		std::cout << " Setting initial field magnitude to: " << sim::H_applied << " T" << std::endl;

		std::cout << " Writing bit sequence ";
      for(auto i=0; i<hamr::internal::bit_sequence.size(); i++){ std::cout << hamr::internal::bit_sequence[i] << " ";}
		std::cout << std::endl;
		std::cout << " Number of bits in x,y:\t" << n_bits_per_tack << "\t" << hamr::internal::num_tracks << std::endl;
		std::cout << " Initial Head position:\t" << hamr::internal::head_position_x*0.1 << "\t" << hamr::internal::head_position_y*0.1 << " nm" << std::endl;
		std::cout << " Head velocity:\t" << hamr::internal::head_speed*1e-10 << "\tm/s" << std::endl;
		std::cout << " Field ramp time:\t" << ramp_time*mp::dt_SI << "\ts" << std::endl;
		std::cout << " Time per bit:\t" << bit_time*mp::dt_SI << "\ts" << std::endl;
		// Check that field ramp time < bit time
		if(2*ramp_time >= bit_time){
         terminaltextcolor(RED);
			std::cerr << "Error - value for \'hamr:field-ramp-time\' too large. Make sure is less than half time per bit." << std::endl;
         terminaltextcolor(WHITE);
			zlog << zTs() << "Error - value for \'hamr:field-ramp-time\' too large. Make sure is less than half time per bit." << std::endl;
         err::vexit();
		}
		std::cout << " Time per track:\t" << track_time*mp::dt_SI << "\ts" << std::endl;
		std::cout << " New total simulated time:\t" << total_time*mp::dt_SI << "\ts" << std::endl;

		int bit_tot = 0;
		int bit = 0;
		int track = 0;
		while(sim::time-sim::equilibration_time < total_time){ // loop over whole time
			uint64_t tmp_time = 0;

			track = int(round(bit_tot/n_bits_per_tack));
			bit = bit_tot - track*n_bits_per_tack;

			hamr::internal::head_position_x = head_position_initial + BL*bit;	
			hamr::internal::head_position_y = TW*(0.5 + track) + padding;
			zlog << zTs() << "New head position:" << hamr::internal::head_position_x*0.1 << ", " << hamr::internal::head_position_y*0.1 << " nm" << std::endl;

			while(tmp_time < bit_time){
				// track = int(bit_tot/n_bits_per_tack);
				// bit = bit_tot - track*n_bits_per_tack;
				hamr::internal::head_position_x = head_position_initial + (speed/**1e-10*/) * tmp_time * mp::dt_SI + BL*bit;
				hamr::internal::head_position_y = TW*(0.5 + track) + padding;
				const double H_app_dir = static_cast<double>(hamr::internal::bit_sequence[bit_tot]);

				// Update applied field value depending on trapezoidal time profile
				const double H_app_abs = fabs(sim::H_applied) + hamr::internal::update_field_time_trapz_profile(tmp_time, ramp_time, bit_time);
				// Determine sign of applied field
				sim::H_applied = H_app_abs * H_app_dir; 

				// Integrate system
				sim::integrate(sim::partial_time);

				// Calculate magnetisation statistics
				stats::update();
				// Output data
				vout::data();

				tmp_time++;
			}
			// Print update on screen and log file
			zlog << zTs() << "Bit " << bit_tot << " (" << bit << "," << track << ") written" << std::endl;

			bit_tot++;
		}

      return;
   } // end of hamr_continuous
   
} // end of hamr namespace