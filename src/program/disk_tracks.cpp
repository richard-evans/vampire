//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Time Series program
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///

// Standard Libraries
#include <iostream>
#include <math.h>

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
#include "cells.hpp"

namespace program{

/// @brief Function to calculate magnetisation over a time series
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
///
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///
namespace track_parameters{

   int num_bits_per_track = 4;
   int num_tracks = 1;
   double fly_height = 100;

   double bit_size = 1000;
   double bit_width = 600;

   double xb = bit_size/2.0;
   double yb = bit_width/2.0;
   double zb = bit_size/2.0;

   double cross_track_velocity = 0.0;
   double down_track_velocity = 0.01;

   double initial_x_position = 0;
   double initial_z_position = -2000;

   double Ms = 1;
   int num_bits = (num_bits_per_track +1)*(num_tracks +1);
   std::vector < double > x_track_array(num_bits,0.0);
   std::vector < double > z_track_array(num_bits,0.0);
   std::vector < double > bit_magnetisation(num_bits,0.0);

   double y_track = -fly_height - bit_width/2.0;

}



void create_tracks(){

   using namespace track_parameters;

   int bit = 0;
   int M = 1.0;
   for (double x = -(num_tracks/2.0)*bit_size; x < (num_tracks/2.0)*bit_size + bit_size; x = x + bit_size){
      for (double z = -(num_bits_per_track/2.0)*bit_size; z < (num_bits_per_track/2.0)*bit_size + bit_size; z = z + bit_size){
         x_track_array[bit] = x;
         z_track_array[bit] = z;
         bit_magnetisation[bit] = M;
         bit++;
         M = M*-1;
      }
      M = M*-1;
   }

}



void calculate_field(int cell,int step){


   using namespace track_parameters;



   double down_track_position = initial_z_position + down_track_velocity*step;
   double cross_track_position = initial_x_position + cross_track_velocity*step;


   sim::track_field_x[cell] = 0.0;
   sim::track_field_y[cell] = 0.0;
   sim::track_field_z[cell] = 0.0;


   //cell position in A
   double x_cell = cross_track_position + cells::pos_and_mom_array[4*cell+0];
   double y_cell = -fly_height          - cells::pos_and_mom_array[4*cell+1];
   double z_cell = down_track_position  + cells::pos_and_mom_array[4*cell+2];
   //std::cout
   //std::cout << x_cell << '\t' << y_cell << '\t' << z_cell << std::endl;
   //loop over all bits to calcualte the field from each bit
   for (int bit = 0; bit < num_bits; bit++){

      std::vector <double > H(3,0.0);
      //bit positions in A
      double x_bit = x_track_array[bit];
      double y_bit = y_track;
      double z_bit = z_track_array[bit];

      //pcalcualtes the prefactor (M/4pi)
      double prefactor = Ms*bit_magnetisation[bit]/(4.0*3.14);
      //std::cout <<prefactor <<std::endl;
      //calcualtes the vector in A from the cell to the bits
      double x = x_cell - x_bit;
      double y = y_cell - y_bit;
      double z = z_cell - z_bit;


      for (int k =1; k < 3; k ++){
         int k1 = pow((-1),k);

         for (int l =1; l < 3; l ++){
            int l1 = pow((-1),l);

            for (int m =1; m < 3; m ++){
               int m1 = pow((-1),m);
               int klm1 = pow((-1),l+k+m);


               double y_l1yb = y+l1*yb;
               double x_k1xb = x+k1*xb;
               double z_m1zb = z+m1*zb;

               if (x_k1xb ==0) x_k1xb = 0.001;
               if (y_l1yb ==0) y_l1yb = 0.001;
               if (z_m1zb ==0) z_m1zb = 0.001;

               double mod = sqrt(y_l1yb*y_l1yb + x_k1xb*x_k1xb);
               double frac = y_l1yb*x_k1xb/mod;
               double sq = x_k1xb*x_k1xb + y_l1yb*y_l1yb + z_m1zb*z_m1zb;

               H[0] += klm1*log(z_m1zb + sqrt(sq));
               H[1] += klm1*frac*atan((x_k1xb*z_m1zb)/(y_l1yb*sqrt(sq)));
               H[2] += klm1*log(x_k1xb + sqrt(sq));

            //   if (cell == 0 ) std::cout << "A" << cell << '\t' << sim::time <<  "\t" << x << '\t' << y << '\t' << z << "\t" << sim::track_field_y[cell]  << std::endl;

            }
         }
      }//std::cin.get();


      sim::track_field_x[cell] += prefactor*H[0];
      sim::track_field_y[cell] += prefactor*H[1];
      sim::track_field_z[cell] += prefactor*H[2];

  //    if (cell == 0 ) std::cout  << sim::time << "\t" <<down_track_position <<   "\t" << sim::track_field_x[cell]  << '\t' << sim::track_field_y[cell]  << '\t' << sim::track_field_z[cell]  << std::endl;
   }

}

void tracks(){

  using namespace track_parameters;


	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::tracks has been called" << std::endl;

	double temp=sim::temperature;

   sim::track_field_x.resize(cells::num_cells);
   sim::track_field_y.resize(cells::num_cells);
   sim::track_field_z.resize(cells::num_cells);





   create_tracks();




	sim::temperature=temp;


  while(sim::time<sim::equilibration_time+sim::total_time){


          for (int lc = 0; lc < cells::num_local_cells; lc++){
             int cell = cells::cell_id_array[lc];
            // std::cout << cell << std::endl;
             calculate_field(cell, sim::time);
          }

    // Integrate system
    sim::integrate(sim::partial_time);

    // Calculate magnetisation statistics
    stats::mag_m();

    // Output data
    vout::data();

       double down_track_position = initial_z_position + down_track_velocity*sim::time;
    std::cout  << sim::time << "\t" <<down_track_position <<   "\t" << sim::track_field_x[0]  << '\t' << sim::track_field_y[0]  << '\t' << sim::track_field_z[0]  << std::endl;

	}

}

}//end of namespace program
