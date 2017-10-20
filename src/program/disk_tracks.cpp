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
#include <fstream>

namespace program{

//------------------------------------------------------------------------------
// Namespace defining dimensions and parameters for tracks
//------------------------------------------------------------------------------

namespace track_parameters{

   // specify number of bits and tracks
   int num_bits_per_track = 1;
   int num_tracks = 1;

   // distance of tracks from read head
   double fly_height = 0.0; // Angstroms

   double bit_size = 1000.0; // size of bits in x-direction (cross track)
   double bit_width = 10000.0; // size of bits in z-direction (down track)
   double bit_depth = 600.0; // depth of bits along y-direction

   double Ms = 0.1;// mu0 Ms in Tesla

   int num_bits; // total number of bits

   std::vector < double > x_track_array; // stores coordinates of bits
   std::vector < double > z_track_array; // stores coordinates of bits
   std::vector < double > bit_magnetisation; // stores magnetization of bits
   //double y_track = -bit_depth*0.5;
}

// typesafe sign function
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


//------------------------------------------------------------------------------
// Function to create tracks
//------------------------------------------------------------------------------

void create_tracks(){
   int M = 1;


   // temporary constants defining half sizes of bits
   const double xb = track_parameters::bit_size*0.5;
   const double yb = track_parameters::bit_depth*0.5;
   const double zb = track_parameters::bit_width*0.5;

   int bit = 0;

   const int start_x = -(track_parameters::num_tracks*xb) + xb;
   const int end_x = (track_parameters::num_tracks*xb) + xb;
   const int bs = track_parameters::bit_size;
   const int bw = track_parameters::bit_width;
   const int start_z = -(track_parameters::num_bits_per_track*zb) + zb;
   const int end_z = (track_parameters::num_bits_per_track*zb) + zb;

    for (int x = start_x; x < end_x; x = x + bs){
       for (double z = -start_z; z < end_z; z = z + bw){
          track_parameters::x_track_array[bit] = x;
          track_parameters::z_track_array[bit] = z;
   //
          track_parameters::bit_magnetisation[bit] = M;
          bit++;
          std::cout << bit << "\t" << track_parameters::num_bits<<std::endl;
   //    //   M = M*-1;
       }
     //  M = M*-1;
    }
}


std::vector <double > calculate_field(double cx, double cy, double cz, int step){

   double down_track_position = sim::initial_down_track_position + sim::down_track_velocity*step;
   double cross_track_position = sim::initial_cross_track_position + sim::cross_track_velocity*step;

   double prefactor = track_parameters::Ms/(4.0*M_PI);

   //cell position in Angstrom
   double x_cell = cross_track_position + cx;
   double y_cell = cy;
   double z_cell = down_track_position  + cz;

   const double xb = track_parameters::bit_size  * 0.5;
   const double yb = track_parameters::bit_depth * 0.5;
   const double zb = track_parameters::bit_width * 0.5;

   const double y_track = -track_parameters::bit_depth/2.0 - track_parameters::fly_height;

   std::vector <double > B(3,0.0);

   for (int bit = 0; bit < track_parameters::num_bits; bit++){

      //bit positions in A
      double x_bit = track_parameters::x_track_array[bit];
      double y_bit = y_track;
      double z_bit = track_parameters::z_track_array[bit];

      //pcalcualtes the prefactor (M/4pi)

      //calculates the vector in A from the cell to the bits
      double x = x_cell - x_bit;
      double y = y_cell - y_bit;
      double z = z_cell - z_bit;


      double Bx = 0.0;
      double By = 0.0;
      double Bz = 0.0;

      for(int k=1; k<4; k++){

          // predefine power as fixed for loop iteration
          const double m1k = pow(-1,k);

          for(int l=1; l<4; l++){

             // predefine power as fixed for loop iteration
             const double m1l = pow(-1,l);

             for(int m=1; m<4; m++){

                const double m1m = pow(-1,m);
                const double m1klm = pow(-1,k+l+m);

                const double xp = x + xb*m1k;
                const double yp = y + yb*m1l;
                const double zp = z + zb*m1m;

                const double r = sqrt( xp*xp + yp*yp + zp*zp );

                const double xabs = fabs(xp);
                const double yabs = fabs(yp);

                Bx += m1klm * log(zp+r);
                By += m1klm * log(xp+r);
                Bz += m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));

             }
          }
      }


      B[0] = Bx*prefactor;
      B[1] = By*prefactor;
      B[2] = Bz*prefactor;
   }

   return B;

}

void tracks(){


     sim::track_field_x.resize(cells::num_cells,0.0);
     sim::track_field_y.resize(cells::num_cells,0.0);
     sim::track_field_z.resize(cells::num_cells,0.0);

   //   for (int cell = 0; cell <cells::num_cells; cell++ ){
   //      sim::track_field_x[cell] = 0.0;
   //      sim::track_field_y[cell] = 0.0;
   //      sim::track_field_z[cell] = 0.0;
   //   }

   // define total number of bits
   track_parameters::num_bits = (track_parameters::num_bits_per_track)*(track_parameters::num_tracks);

   track_parameters::x_track_array.resize(track_parameters::num_bits,0.0);
   track_parameters::z_track_array.resize(track_parameters::num_bits,0.0);
   track_parameters::bit_magnetisation.resize(track_parameters::num_bits,0.0);

   create_tracks();


std::vector <double > B(3,0.0);

  while(sim::time<sim::equilibration_time+sim::total_time){

     for (int lc = 0; lc < cells::num_local_cells; lc++){
        int cell = cells::cell_id_array[lc];

          const double cx = cells::pos_and_mom_array[4*cell+0];
          const double cy = cells::pos_and_mom_array[4*cell+1];
          const double cz = cells::pos_and_mom_array[4*cell+2];
          B = calculate_field(cx, cy, cz, sim::time);
          sim::track_field_x[cell] = B[0];
          sim::track_field_y[cell] = B[1];
          sim::track_field_z[cell] = B[2];
       }


    // Integrate system
    sim::integrate(sim::partial_time);

    // Calculate magnetisation statistics
    stats::mag_m();

    // Output data
    vout::data();

	}

}
}//end of namespace program
