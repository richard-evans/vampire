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

		std::cout << " >>>> Performing HAMR continuous simulation" << std::endl;
		zlog << zTs() << " >>>> Performing HAMR continuous simulation" << std::endl;
		
		// Calibrate head region to be not larger than system size
		if(hamr::internal::H_bounds_x > hamr::internal::system_dimensions[0]){
			hamr::internal::H_bounds_x = hamr::internal::system_dimensions[0];
			hamr::internal::bit_spacing_x = 0.0;
			std::cout << " >>>> Resizing sim:hamr-head-field-x to\t" << hamr::internal::H_bounds_x << "\tAng" << std::endl;
			zlog << zTs() << " >>>> Resizing sim:hamr-head-field-x to\t" << hamr::internal::H_bounds_x << "\tAng" << std::endl;
		}
		if(hamr::internal::H_bounds_y > hamr::internal::system_dimensions[1]){
			hamr::internal::H_bounds_y = hamr::internal::system_dimensions[1];
			hamr::internal::bit_spacing_y = 0.0;
			std::cout << " >>>> Resizing sim:hamr-head-field-y to\t" << hamr::internal::H_bounds_y << "\tAng" << std::endl;
			zlog << zTs() << " >>>> Resizing sim:hamr-head-field-y to\t" << hamr::internal::H_bounds_y << "\tAng" << std::endl;
		}

		// Initialise number of bits so that at least one bit is written
		int N_bit_x = 1;
		int N_bit_y = 1;
		int x_bit = 0;
		int y_bit = 0;
		int bit = 0;
		// Initial position is with head away from film
		std::cout << "H_bounds_x\t" << hamr::internal::H_bounds_x << std::endl;
		std::cout << "H_bounds_y\t" << hamr::internal::H_bounds_y << std::endl;
		hamr::internal::head_position[0] = 0.0 - hamr::internal::H_bounds_x*0.5;
		hamr::internal::bit_spacing_x = 0.0;
		hamr::internal::head_position[1] = hamr::internal::H_bounds_y*0.5 + hamr::internal::bit_spacing_y;

		// If simulation of single bit, set the head in the centre of the y-system dimension
		if(hamr::internal::single_bit==true){
			std::cout << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
			zlog << zTs() << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
			hamr::internal::H_osc_amplit = hamr::internal::H_bounds_x;
			std::cout << " >>>> Setting hamr:field-oscillation-frequency to\t" << hamr::internal::H_osc_amplit << "\tAng" << std::endl;
			zlog << zTs() << " >>>> Setting hamr:field-oscillation-frequency to\t" << hamr::internal::H_osc_amplit << "\tAng" << std::endl;
//			hamr::internal::head_position[1] = hamr::internal::system_dimensions[1]*0.5;
		}
		// Otherwise, determine the total number of bits in x and y-directions
		else{
			// Determine number of bits
			if( floor(hamr::internal::system_dimensions[0]/(hamr::internal::H_bounds_x+hamr::internal::bit_spacing_x + 1.0))>0 ){
				N_bit_x = floor(hamr::internal::system_dimensions[0]/(hamr::internal::H_bounds_x+hamr::internal::bit_spacing_x + 1.0)); // number of bits in x
			}
			if( floor(hamr::internal::system_dimensions[1]/(hamr::internal::H_bounds_y+hamr::internal::bit_spacing_y + 1.0))>0 ){
				N_bit_y = floor(hamr::internal::system_dimensions[1]/(hamr::internal::H_bounds_y+hamr::internal::bit_spacing_y + 1.0)); // number of bits in y
			}
		}

		// Set times in terms of integragion steps
		const uint64_t peak_time = int( (((hamr::internal::H_bounds_x + hamr::internal::bit_spacing_x*0.5)/(hamr::internal::head_speed))/6.0) /mp::dt_SI );
		const uint64_t write_time = 6*peak_time;
		const uint64_t pre_write_time = 0.5*write_time;
		const uint64_t ramp_time = int(hamr::internal::H_ramp_time/mp::dt_SI);
		std::cout << " head_speed\t" << hamr::internal::head_speed*1e-10 << "\tm/s" << std::endl;
		std::cout << " peak_time\t"  << peak_time*mp::dt_SI   << "\ts  " << std::endl;
		std::cout << " write_time\t" << write_time*mp::dt_SI  << "\ts  " << std::endl;
		std::cout << " ramp_time\t"  << ramp_time*mp::dt_SI   << "\ts  " << std::endl;
		std::cout << " equl time\t"  << sim::equilibration_time*mp::dt_SI << "\ts" << std::endl;
		std::cout << " total time\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time )*mp::dt_SI << "\ts" << std::endl;
		std::cout << " total time + equl time\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;

		const double Hmax = hamr::internal::Hmax; // max field
		const double Hmin = hamr::internal::Hmin;
		sim::H_applied = Hmin;

		std::cout << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
		std::cout << " >>>> Initial Head position:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << "\tAng" << std::endl;
		std::cout << " >>>> Head velocity:\t" << hamr::internal::head_speed*1e-10 << "\tm/s" << std::endl;
		std::cout << " >>>> Time per bit:\t" << (write_time*N_bit_x*N_bit_y)*mp::dt_SI << "\ts" << std::endl;
		std::cout << " >>>> New total simulated time:\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;
		zlog << zTs() << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
		zlog << zTs() << " >>>> Initial Head poistion:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << "\tAng" << std::endl;
		zlog << zTs() << " >>>> Head velocity:\t" << hamr::internal::head_speed*1e-10 << "\tm/s" << std::endl;
		zlog << zTs() << " >>>> Time per bit:\t" << (write_time*N_bit_x*N_bit_y)*mp::dt_SI << "\ts" << std::endl;
		zlog << zTs() << " >>>> New total simulated time:\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;

//		while(sim::time-sim::equilibration_time<write_time*N_bit_x*N_bit_y && y_bit<N_bit_y && hamr::internal::head_position[1]<hamr::internal::system_dimensions[1]){
		while(sim::time-sim::equilibration_time < write_time*N_bit_x*N_bit_y+pre_write_time*2 &&
				hamr::internal::head_position[1]*y_bit<hamr::internal::system_dimensions[1]){  /*&& bit<N_bit_x*N_bit_y+1 && y_bit<N_bit_y*/
			// Update head position in y-direction
			hamr::internal::head_position[1] = hamr::internal::H_bounds_y*(0.5+y_bit) + hamr::internal::bit_spacing_y;

//			while( hamr::internal::head_position[0]-0.5*hamr::internal::H_bounds_x < hamr::internal::system_dimensions[0] &&
//						 hamr::internal::head_position[0]-0.5*hamr::internal::H_bounds_x < N_bit_x*(hamr::internal::H_bounds_x+hamr::internal::bit_spacing_x) &&
//						 bit<(N_bit_x)*N_bit_y+1 ){
			while( hamr::internal::head_position[0] < hamr::internal::system_dimensions[0]+0.5*hamr::internal::H_bounds_x ){ // &&
//					hamr::internal::head_position[0]  -0.5*hamr::internal::H_bounds_x   < N_bit_x*(hamr::internal::H_bounds_x+hamr::internal::bit_spacing_x) ){

				std::cout << "\n >>>> Moving HEAD **** New HEAD position:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << "\tAng\n" << std::endl;
				zlog << zTs() << " >>>> Moving HEAD **** New HEAD position:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << "\tAng" << std::endl;

//				while( hamr::internal::head_position[0] < hamr::internal::system_dimensions[0]+0.5*hamr::internal::H_bounds_x &&
//						 hamr::internal::head_position[0] < hamr::internal::H_bounds_x*(bit - N_bit_x*(y_bit)+1)    &&
//						 x_bit < N_bit_x ){
				while( hamr::internal::head_position[0] < hamr::internal::H_bounds_x*(x_bit+1) &&
						 hamr::internal::head_position[0] < hamr::internal::system_dimensions[0]+0.5*hamr::internal::H_bounds_x){
					// Update head position in x-direction
					hamr::internal::head_position[0] = ( (hamr::internal::head_speed) * (sim::time-sim::equilibration_time)*mp::dt_SI - 0.5*hamr::internal::H_bounds_x ) - (N_bit_x*hamr::internal::H_bounds_x)*y_bit;

					// loop over partial time
					for(int tt=0; tt < sim::partial_time; tt++){
						const uint64_t field_time = (sim::time-sim::equilibration_time) - write_time*bit;
						const uint64_t ramp_time_end  = write_time - ramp_time;
						// Update applied field value depending on trapezoidal time profile
      				hamr::internal::update_field_time_trapz_profile(field_time, ramp_time, ramp_time_end, pre_write_time, write_time, sim::H_applied);
						// Integrate system
						sim::integrate(1);
					} // end loop over time-step-increment=partial_time
					// Calculate magnetisation statistics
					stats::mag_m();
					// Output data
					vout::data();
				}  // End of one single bit
				std::cout << "\n >>>> Finished writing bit\t(" << x_bit << "," << y_bit << ")\t Hz\t"
						    << sim::H_vec[2]*((-1.0)*double(2*(int(hamr::internal::head_position[0]/hamr::internal::H_osc_amplit)%2)-1)) << "\n" << std::endl;
				zlog << zTs() << " >>>> Finished writing bit\t(" << x_bit << "," << y_bit << ")\t Hz\t"
						        << sim::H_vec[2]*((-1.0)*double(2*(int(hamr::internal::head_position[0]/hamr::internal::H_osc_amplit)%2)-1)) << std::endl;
				++x_bit;
				++bit;
			}  // End of continuous displacement along x
			std::cout << " >>>> Final track head position:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << std::endl;
			std::cout << " \t     N_bit_x (hamr::internal::H_bounds_x + hamr::internal::bit_spacing_x)\t" << N_bit_x*(hamr::internal::H_bounds_x+hamr::internal::bit_spacing_x) << std::endl;
			std::cout << " \t     x_bit\t" << x_bit << std::endl;
			std::cout << " \t     y_bit\t" << y_bit << std::endl;
			std::cout << " \t     bit\t" << bit << std::endl;
			std::cout << " \t     N_bit_x\t" << N_bit_x << std::endl;
			std::cout << " \t     N_bit_y\t" << N_bit_y << std::endl;
			std::cout << "\n >>>> Reset head at the beginning of downtrack and move along offtrack direction\n" << std::endl;
			zlog << zTs() << " >>>> Final track head position:\t" << hamr::internal::head_position[0] << "\t" << hamr::internal::head_position[1] << std::endl;
			zlog << zTs() << " >>>> Reset head at the beginning of downtrack and move along offtrack direction" << std::endl;
			// Update head position to allow correct evaluation of statement in while() loops
			hamr::internal::head_position[0] = 0.0 - 0.5*hamr::internal::H_bounds_x;
			hamr::internal::head_position[1] += hamr::internal::H_bounds_y + hamr::internal::bit_spacing_y;
			// Reset counter for bit in downtrack direction
			x_bit = 0;
			++y_bit;
		}  // End of loop over N_bit_y

//			// force outputting last point of simulation
//			if(sim::time-sim::equilibration_time>=write_time*N_bit_x*N_bit_y){
//				// Disable laser
//				sim::head_laser_on=false;
//				std::cout << "\n>>>> Disable laser and integrate system for 1 time-step" << std::endl;
//				// Integrate
//				sim::integrate(1);
//				// Calculate magnetisation statistics
//				stats::mag_m();
//				// Output data
//				vout::data();
//				std::cout << "\n>>>> Outputting system at the end of HAMR continuous simulations\n" << std::endl;
//				// Reactivate laser
//				sim::head_laser_on=true;
//			}

      return;
   } // end of hamr_continuous
   
} // end of hamr namespace