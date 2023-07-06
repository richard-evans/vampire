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

		//const double Tmin = hamr::internal::Tmin;
		// Determine times in terms of integration steps
		const double BL = hamr::internal::bit_size;
		const double TW = hamr::internal::track_size;
		const double speed = hamr::internal::head_speed;
		const double NPS = hamr::internal::NPS;
		const double padding = hamr::internal::track_padding;
		//const int n_bits = hamr::internal::num_bits;
		const int n_bits_per_track = hamr::internal::bits_per_track;
		const int n_tracks = hamr::internal::num_tracks;
		const uint64_t NPS_time   = int(round((NPS/speed)/mp::dt_SI)); // time to cover extra sweeping due to NPS
		const uint64_t bit_time   = int(round((BL/speed)/mp::dt_SI));
		const uint64_t extra_time = int(round(((BL*0.5)/speed)/mp::dt_SI)) + NPS_time;  // time to sweep across half of bit size necessary to cover initial distance of head from track
		const uint64_t track_time = n_bits_per_track*bit_time;
		const uint64_t total_track_time = track_time + extra_time*2;
		const uint64_t total_time = n_tracks * total_track_time;
		const uint64_t final_time = int(1.0e-11/mp::dt_SI);  // Integrate for an extra ps to allow all regions to equilibrate to Tmin
		const uint64_t rise_time = int(round(hamr::internal::H_rise_time/mp::dt_SI));
		const uint64_t fall_time = int(round(hamr::internal::H_fall_time/mp::dt_SI));
		const double Deltax = speed * mp::dt_SI; // speed is units of Angstrom/second
		// Initial head position
		const double head_position_initial = -BL*0.5 - NPS;
		hamr::internal::head_position_x = head_position_initial;
		hamr::internal::head_position_y = 0.5*hamr::internal::system_dimensions_y;
		// Initialise field magnitude to min value at beginning of simulation
		sim::H_applied = hamr::internal::Hmin;
		std::cout << " Setting initial field magnitude to: " << sim::H_applied << " T" << std::endl;

		std::cout << " Bit sequence to be written: ";
		for(size_t i=0; i<hamr::internal::bit_sequence.size(); i++){ std::cout << hamr::internal::bit_sequence[i] << " ";}
		std::cout << std::endl;
		std::cout << " Bits to be written: " << n_bits_per_track*n_tracks << ", in x,y: " << n_bits_per_track << ", " << n_tracks << std::endl;
		zlog << zTs() << "Bits to be written: " << n_bits_per_track*n_tracks << ", in x,y: " << n_bits_per_track << ", " << n_tracks << std::endl;

		// Print to screen hamr writing setup parameters
		std::cout << " Head velocity: " << hamr::internal::head_speed*1e-10 << " m/s" << std::endl;
		std::cout << " Time per bit: " << bit_time*mp::dt_SI << " s" << std::endl;
		// Check that field ramp time < bit time
		if(rise_time+fall_time >= bit_time){
         terminaltextcolor(RED);
			std::cerr << "Error - value for \'hamr:field-rise-time\' " << rise_time*mp::dt_SI << " s or \'hamr:field-fall-time\' " << fall_time*mp::dt_SI << " s is too large. Make sure is less than half time per bit: " << bit_time*mp::dt_SI << " s." << std::endl;
         terminaltextcolor(WHITE);
			zlog << zTs() << "Error - value for \'hamr:field-rise-time\' " << rise_time*mp::dt_SI << " s or \'hamr:field-fall-time\' " << fall_time*mp::dt_SI << " s is too large. Make sure is less than half time per bit: " << bit_time*mp::dt_SI << " s." << std::endl;
         err::vexit();
		}
		std::cout << " Time per track: " << total_track_time*mp::dt_SI << " s" << std::endl;
		std::cout << " New total simulation time: " << (total_time + final_time)*mp::dt_SI << " s" << std::endl;

		//--------------------------------------------//
		// Start hamr simulation
		//--------------------------------------------//
		int bit_tot = 0;
		int bit = 0;
		int track = 0;
		const uint64_t total_sim_time = total_time + sim::equilibration_time;
		const uint64_t total_final_time = final_time + total_sim_time;
		sim::total_time = total_final_time - sim::equilibration_time; // update global value of total time, needed for output of atomic configurations
		while(sim::time < total_sim_time){ // loop over whole time

			track = int(round(bit_tot/n_bits_per_track));
			bit = bit_tot - track*n_bits_per_track;
			// Set head position at beginning of track
			hamr::internal::head_position_x +=Deltax; //= head_position_initial + BL*bit;
			hamr::internal::head_position_y = TW*(0.5 + track) + padding;
			zlog << zTs() << "Head position at beginning of track " << track+1 << " : " << hamr::internal::head_position_x*0.1 << " nm, " << hamr::internal::head_position_y*0.1 << " nm" << std::endl;

			//-----------------------------------------//
			// Head starts half bit away from track
			// magnetic field off
			//-----------------------------------------//
			// Decleare integer time step within track
			uint64_t tmp_track_time = 0;
			while(tmp_track_time < extra_time){

				// Switch off external field
				sim::H_applied = 0.0;

				// Update head position in downtrack
				hamr::internal::head_position_x += Deltax;

				// Integrate system
				sim::integrate(sim::partial_time);
				// Calculate magnetisation statistics
				stats::update();
				// Output data
				vout::data();

				tmp_track_time++;
			} // end of loop over head outside track region
			zlog << zTs() << "Head position at beginning of first bit in track " << track+1 << " : " << (hamr::internal::head_position_x+Deltax)*0.1 << " nm, " << hamr::internal::head_position_y*0.1 << " nm" << std::endl;

			//-----------------------------------------//
			// Enter in track => within bits
			//-----------------------------------------//
			while(tmp_track_time < track_time){

				// Decleare integer time step within bit
				uint64_t tmp_bit_time = 0;
				while(tmp_bit_time < bit_time){

					// Move head along downtrack
					hamr::internal::head_position_x += Deltax;

					// Determine field polarisation within bit
					const double H_app_dir = static_cast<double>(hamr::internal::bit_sequence[bit_tot]);
					// Update applied field value depending on trapezoidal time profile
					const double H_app_abs = fabs(sim::H_applied) + hamr::internal::update_field_time_trapz_profile(tmp_bit_time, rise_time, fall_time, bit_time);
					// Determine sign of applied field
					sim::H_applied = H_app_abs * H_app_dir;

					// Integrate system
					sim::integrate(sim::partial_time);

					// Calculate magnetisation statistics
					stats::update();
					// Output data
					vout::data();

					tmp_bit_time++;
					tmp_track_time++;
				} // end of loop over single bit region

				// Print update on screen and log file
				zlog << zTs() << "Written bit " << bit_tot+1 << " : bit " << bit+1 << " on track " << track+1 <<  std::endl;
				zlog << zTs() << "Head moving to position: " << (hamr::internal::head_position_x+Deltax)*0.1 << " nm, " << hamr::internal::head_position_y*0.1 << " nm" << std::endl;

				bit_tot++;
			} // end of loop over track region

			//--------------------------------------------------//
			// Let head go beyond track edge for half bit size,
			// with magnetic field off, to ensure
			// last bit of track is well written
			//--------------------------------------------------//
			while(tmp_track_time < total_track_time){

				// Switch off external field
				sim::H_applied = 0.0;

				// Update head position in downtrack
				hamr::internal::head_position_x += Deltax;

				// Integrate system
				sim::integrate(sim::partial_time);
				// Calculate magnetisation statistics
				stats::update();
				// Output data
				vout::data();

				tmp_track_time++;
			} // end of loop over head outside track region
			zlog << zTs() << "Head position at end of track " << track+1 << " : " << (hamr::internal::head_position_x+Deltax)*0.1 << ", " << hamr::internal::head_position_y*0.1 << " nm" << std::endl;

			// Reset head position in downtrack once a track has been written
			hamr::internal::head_position_x = head_position_initial;

		}

		zlog << zTs() << "Disabling laser and external field and integrating system for extra 10 ps" << std::endl;
		while(sim::time < total_final_time){
			// Disable laser
			hamr::head_laser_on=false;
			// Switch off external field
			sim::H_applied = 0.0;
			// Set system temperature as minimum temperature
			sim::temperature=sim::Tmin;

			// Integrate system
			sim::integrate(sim::partial_time);
			// Calculate magnetisation statistics
			stats::update();
			// Output data
			vout::data();
		}

      return;
   } // end of hamr_continuous

} // end of hamr namespace
