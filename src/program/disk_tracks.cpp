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
#include <fstream>
#include <sstream>
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
#include "micromagnetic.hpp"
#include "../environment/internal.hpp"
#include <fstream>


namespace program{

//------------------------------------------------------------------------------
// Namespace defining dimensions and parameters for tracks
//------------------------------------------------------------------------------



namespace tp{

   // specify number of bits and tracks
   int num_bits_per_track;
   int num_tracks;



   // distance of tracks from read head
   double fly_height; // Angstroms

   double bit_size; // size of bits in x-direction (cross track)
   double bit_width; // size of bits in z-direction (down track)
   double bit_depth; // depth of bits along y-direction

   double track_gap;
   double bit_gap;

   double total_bit_size;
   double total_bit_width;
   double total_bit_depth;

   double Ms;// mu0 Ms in Tesla

   int num_bits; // total number of bits

   std::vector < double > x_track_pos; // stores coordinates of bits
   std::vector < double > z_track_pos; // stores coordinates of bits
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

   std::ifstream ifile;
   ifile.open("track_ms");

   int i = 0;

   int track_num;
   int bit_num;

   std::vector <double > bitms(tp::num_tracks*tp::num_bits_per_track,0.0);

   if (sim::track_ms_file == true){

     // read number of atoms in this file
     std::string line; // declare a string to hold line of text
     while(getline(ifile,line) ){
       std::stringstream line_stream(line);
       line_stream >> track_num >> bit_num >> M;
       std::cout << (track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << '\t' << bitms.size() << "\t" << tp::num_bits_per_track << '\t' <<tp::num_tracks<<std::endl;
       bitms[(track_num-1)*(tp::num_bits_per_track) + (bit_num-1)] = M*tp::Ms;
       std::cout <<(track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << std::endl;
     }
   }
   else {
  //   std::cout << "ENTER2" << std::endl;
     int M = 1;
     for (int bit  =0; bit < tp::num_bits; bit++){
       for (int track =0; track < tp::num_tracks; track++){
         bitms[(track-1)*(tp::num_bits_per_track) + (bit-1)] = M;
         M *=-1;
       }
     }
   }
   std::cout << "A" << std::endl;

   // temporary constants defining half sizes of bits
   const double xb = tp::bit_size*0.5;
   const double yb = tp::bit_depth*0.5;
   const double zb = tp::bit_width*0.5;

   int bit = 0;

   const int start_x = -(tp::num_tracks*xb) + xb;
   const int end_x = (tp::num_tracks*xb) + xb;
   const int bs = tp::total_bit_size;
   const int bw = tp::total_bit_width;
   const int start_z = -(tp::num_bits_per_track*zb) + zb;
   const int end_z = (tp::num_bits_per_track*zb) + zb;

    for (int x = start_x; x < end_x; x = x + bs){
       for (double z = start_z; z < end_z; z = z + bw){
          tp::x_track_pos[bit] = x;
          tp::z_track_pos[bit] = z;
          tp::bit_magnetisation[bit] = bitms[bit];
          std::cout << x << '\t' << z <<  "\t" << bw << '\t' << bs << std::endl;
      //    std::cout << bitms[bit] <<std::endl;
    //      std::cout << bit << '\t' <<  tp::bit_magnetisation[bit] <<std::endl;
          bit++;

//      M = M*-1;
       }
  //     M = M*-1;
    }
    std::cout << "A" << std::endl;


//    std::cin.get();
}


std::vector <double > calculate_field(double cx, double cy, double cz, int step){

   double down_track_position = sim::initial_down_track_position + sim::down_track_velocity*step;
   double cross_track_position = sim::initial_cross_track_position + sim::cross_track_velocity*step;

   double prefactor = tp::Ms/(4.0*M_PI);

   //cell position in Angstrom
   double x_cell = cross_track_position + cx;
   double y_cell = cy;
   double z_cell = down_track_position  + cz;

   const double xb = tp::bit_size  * 0.5;
   const double yb = tp::bit_depth * 0.5;
   const double zb = tp::bit_width * 0.5;

   const double y_track = -tp::bit_depth/2.0 - tp::fly_height;

   std::vector <double > B(3,0.0);

   for (int bit = 0; bit < tp::num_bits; bit++){

      //bit positions in A
      double x_bit = tp::x_track_pos[bit];
      double y_bit = y_track;
      double z_bit = tp::z_track_pos[bit];


      //calculates the vector in A from the cell to the bits
      double x = sqrt((x_cell - x_bit)*(x_cell - x_bit));
      double y = sqrt((y_cell - y_bit)*(y_cell - y_bit));
      double z = sqrt((z_cell - z_bit)*(z_cell - z_bit));


    //  std::cout << "distance" << "\t" << x << "\t" << y <<  "\t'" << z << "\t" << tp::bit_magnetisation[bit] <<std::endl;
      double Bx = 0.0;
      double By = 0.0;
      double Bz = 0.0;

      for(int k=1; k<3; k++){

          // predefine power as fixed for loop iteration
          const double m1k = pow(-1,k);

          for(int l=1; l<3; l++){

             // predefine power as fixed for loop iteration
             const double m1l = pow(-1,l);

             for(int m=1; m<3; m++){

                const double m1m = pow(-1,m);
                const double m1klm = pow(-1,k+l+m);

                const double xp = x + xb*m1k;
                const double yp = y + yb*m1l;
                const double zp = z + zb*m1m;
               // std::cout
                const double xabs = fabs(xp);
                const double yabs = fabs(yp);

                double r = sqrt(xp*xp + yp*yp + zp*zp);

                Bx = Bx + m1klm* log(zp + r);
                By = By + m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));
                Bz = Bz + m1klm* log(xp + r);
      //          std::cout << bit  << '\t' << Bx << '\t' << By << '\t' << Bz << std::endl;


             }
          }
      }
  // std::cout <<"field1:\t" <<  Bx << '\t' <<  By <<std::endl;
  // std::cout << "\t" << std::endl;
      B[0] = B[0] + Bx*prefactor*tp::bit_magnetisation[bit];
      B[1] = B[1] + By*prefactor*tp::bit_magnetisation[bit];
      B[2] = B[2] + Bz*prefactor*tp::bit_magnetisation[bit];
    //  std::cout <<"\t\t\t\t\tfield:\t" << Bx << '\t' << By << '\t' << Bz << '\t' << prefactor << '\t' << tp::bit_magnetisation[bit]  << '\t' <<  B[0] << '\t' <<  B[1] << '\t' << B[2] <<std::endl;

   }

   return B;

}

void tracks(){


     tp::num_bits_per_track = sim::track_num_bits_per_track;
     tp::num_tracks = sim::track_num_tracks;



     // distance of tracks from read head
     tp::fly_height = sim::track_fly_height; // Angstroms

     tp::bit_size =  sim::track_bit_size; // size of bits in x-direction (cross track)
     tp::bit_width = sim::track_bit_width; // size of bits in z-direction (down track)
     tp::bit_depth = sim::track_bit_depth; // depth of bits along y-direction

     tp::track_gap = sim::track_track_gap;
     tp::bit_gap = sim::track_bit_gap;

     tp::total_bit_size = tp::bit_size + tp::bit_gap;
     tp::total_bit_width = tp::bit_width + tp::track_gap;
     tp::total_bit_depth = tp::bit_depth;

     tp::Ms = sim::track_Ms;// mu0 Ms in Tesla

     sim::track_field_x.resize(cells::num_cells,0.0);
     sim::track_field_y.resize(cells::num_cells,0.0);
     sim::track_field_z.resize(cells::num_cells,0.0);

   // define total number of bits
    tp::num_bits = (tp::num_bits_per_track)*(tp::num_tracks);


    tp::x_track_pos.resize(tp::num_bits,0.0);
    tp::z_track_pos.resize(tp::num_bits,0.0);
    tp::bit_magnetisation.resize(tp::num_bits,0.0);

    std::vector <double > B(3,0.0);

    int M = 1;

    std::ifstream ifile;
    ifile.open("track_ms");

    int i = 0;

    int track_num;
    int bit_num;
    int bits = 0;
    std::vector <double > bitms(tp::num_tracks*tp::num_bits_per_track,0.0);

    if (sim::track_ms_file == true){
      std::cout << "Creating track from file" << std::endl;

     // read number of atoms in this file
     std::string line; // declare a string to hold line of text
     while(getline(ifile,line) ){
       std::stringstream line_stream(line);
       bits++;
       if (bits > tp::num_bits) {
         std::cout << "error number of bits in file more than number of specified bits" << std::endl;
         err::vexit();
       }
       line_stream >> track_num >> bit_num >> M;
       std::cout << (track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << '\t' << bitms.size() << "\t" << tp::num_bits_per_track << '\t' <<tp::num_tracks<<std::endl;
       if (track_num-1 > tp::num_tracks || bit_num-1  > tp::num_bits_per_track){
         std::cout << "error number of bits in file more than number of specified bits" << std::endl;
         err::vexit();
       }
       bitms[(track_num-1)*(tp::num_bits_per_track) + (bit_num-1)] = M*tp::Ms;
       std::cout <<(track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << std::endl;
     }
   }
   else {
           std::cout << "Creating track using random Ms values" << std::endl;
  //   std::cout << "ENTER2" << std::endl;
     int M = 1;
     for (int bit  =0; bit < tp::num_bits; bit++){
       for (int track =0; track < tp::num_tracks; track++){
         bitms[(track-1)*(tp::num_bits_per_track) + (bit-1)] = M;
         M *=-1;
       }
     }
   }

   // temporary constants defining half sizes of bits
   const double xb = tp::bit_size*0.5;
   const double yb = tp::bit_depth*0.5;
   const double zb = tp::bit_width*0.5;

   int bit = 0;

   const int start_x = -(tp::num_tracks*xb) + xb;
   const int end_x = (tp::num_tracks*xb) + xb;
   const int bs = tp::total_bit_size;
   const int bw = tp::total_bit_width;
   const int start_z = -(tp::num_bits_per_track*zb) + zb;
   const int end_z = (tp::num_bits_per_track*zb) + zb;

    for (int x = start_x; x < end_x; x = x + bs){
       for (double z = start_z; z < end_z; z = z + bw){
          tp::x_track_pos[bit] = x;
          tp::z_track_pos[bit] = z;
          tp::bit_magnetisation[bit] = bitms[bit];
          std::cout << x << '\t' << z <<  "\t" << bw << '\t' << bs <<  "\t" <<bitms[bit] << std::endl;
      //    std::cout << bitms[bit] <<std::endl;
    //      std::cout << bit << '\t' <<  tp::bit_magnetisation[bit] <<std::endl;
          bit++;

       }
    }

   if (environment::internal::LFA_scan == false){

   for (int lc = 0; lc < cells::num_local_cells; lc++){
    int cell = cells::cell_id_array[lc];

    const double cx = cells::pos_and_mom_array[4*cell+0];
    const double cy = cells::pos_and_mom_array[4*cell+1];
    const double cz = cells::pos_and_mom_array[4*cell+2];
    B = calculate_field(cx,cy,cz,sim::time);

    sim::track_field_x[cell] = B[0];
    sim::track_field_y[cell] = B[1];
    sim::track_field_z[cell] = B[2];

  }

   std::ofstream ofile;
   ofile.open ("position.txt");
   int step = 0;

   while(sim::time <sim::equilibration_time+sim::total_time){

   double cross_track_position = sim::initial_cross_track_position + sim::cross_track_velocity*step;
   double down_track_position = sim::initial_down_track_position + sim::down_track_velocity*step;


   for (int lc = 0; lc < cells::num_local_cells; lc++){
     int cell = cells::cell_id_array[lc];

      const double cx = cells::pos_and_mom_array[4*cell+0];
      const double cy = cells::pos_and_mom_array[4*cell+1];
      const double cz = cells::pos_and_mom_array[4*cell+2];

      B = calculate_field(cx,cy,cz,step);
      sim::track_field_x[cell] = B[0];
      sim::track_field_y[cell] = B[1];
      sim::track_field_z[cell] = B[2];
          //   ofile << cell << '\t' << sim::time << '\t' << down_track_position << "\t" << cross_track_position << '\t' << micromagnetic::MR_resistance << "\t" << B[0]<< '\t' << B[1] << '\t' << B[2] <<  std::endl;
         //    std::cout  << sim::time << '\t' << down_track_position << "\t" << cross_track_position << '\t' << micromagnetic::MR_resistance << "\t" << B[0] << '\t' << B[1] << '\t' << B[2] <<  std::endl;
    }

         // Integrate system
  sim::integrate(sim::partial_time);

  ofile << sim::time << '\t' << down_track_position << "\t" << cross_track_position << '\t' << micromagnetic::MR_resistance << "\t" << B[0] << '\t' << B[1] << '\t' << B[2] <<  std::endl;
  step++;

    // Calculate magnetisation statistics
    stats::mag_m();

    // Output data
    vout::data();


	}
 }
//
// else {
//
//
//   while(tp::Ms > -0.11){
//
//
//
//    for (int lc = 0; lc < cells::num_local_cells; lc++){
//      int cell = cells::cell_id_array[lc];
//
//         const double cx = cells::pos_and_mom_array[4*cell+0];
//         const double cy = cells::pos_and_mom_array[4*cell+1];
//         const double cz = cells::pos_and_mom_array[4*cell+2];
//
//             B = calculate_field(cx,cy,cz,0);
//              sim::track_field_x[cell] = B[0];
//              sim::track_field_y[cell] = B[1];
//              sim::track_field_z[cell] = B[2];
//
//          }
//
// for (int i = 0; i < 500; i++){
//     // Integrate system
//     sim::integrate(1000);
//
//     // Calculate magnetisation statistics
//     stats::mag_m();
//
//     // Output data
//     vout::data();
// }
//    //ofile << sim::time << '\t' << tp::Ms << "\t" << micromagnetic::MR_resistance << "\t" << B[0] << '\t' << B[1] << '\t' << B[2] <<  std::endl;
//
//     tp::Ms = tp::Ms - 0.01;
//
// 	}
//
//
//    while(tp::Ms < 0.11){
//
//
//
//     for (int lc = 0; lc < cells::num_local_cells; lc++){
//       int cell = cells::cell_id_array[lc];
//
//          const double cx = cells::pos_and_mom_array[4*cell+0];
//          const double cy = cells::pos_and_mom_array[4*cell+1];
//          const double cz = cells::pos_and_mom_array[4*cell+2];
//
//              B = calculate_field(cx,cy,cz,0);
//               sim::track_field_x[cell] = B[0];
//               sim::track_field_y[cell] = B[1];
//               sim::track_field_z[cell] = B[2];
//
//           }
//
//  for (int i = 0; i < 500; i++){
//      // Integrate system
//      sim::integrate(1000);
//
//      // Calculate magnetisation statistics
//      stats::mag_m();
//
//      // Output data
//      vout::data();
    //}
//
//
//   std::cout << sim::time << '\t' << tp::Ms << "\t" << micromagnetic::MR_resistance << "\t" << B[0] << '\t' << B[1] << '\t' << B[2] <<  std::endl;
//
//      tp::Ms = tp::Ms + 0.01;
//
// 	}
//}
}
}//end of namespace program
