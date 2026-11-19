//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Time Series program
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///   Created:    30/03/2011
///   Revision:   --
///=====================================================================================
///

// Standard Libraries
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
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


   int num_bits; // total number of bits

   std::vector < double > x_track_pos; // stores coordinates of bits
   std::vector < double > z_track_pos; // stores coordinates of bits
   std::vector < double > bit_magnetisation; // stores magnetization of bits
   //double y_track = -bit_depth*0.5;
}

// typesafe sign function
// template <typename T> int sign(T val) {
//     return (T(0) < val) - (val < T(0));
// }

double sign(double a){
      if(a<0.0) return -1.0;
      else return 1.0;
}


//------------------------------------------------------------------------------
// Function to calcualte fields
//------------------------------------------------------------------------------



std::vector <double > calculate_field(double cx, double cy, double cz, int step){

   program::track_pos_z = program::initial_down_track_position  + program::down_track_velocity*step;
   program::track_pos_x = program::initial_cross_track_position + program::cross_track_velocity*step;

   //std::cout << program::track_pos_x << "\t" << program::cross_track_velocity << '\t' << step << std::endl;

   //prefactor includes a conversion from Telsa  to A/m (dived by 4*pi*1e-7) and a divide by 4 pi from the equation in the prism paper
   //this converts the output field to Oe - oh noooo
   double prefactor = program::Ms/(4*M_PI);//program::Ms;//(4.0*M_PI * 4*M_PI* 1e-7);
   //convert Oe to Tesla by diving by 1e4
   //prefactor = prefactor*1e-4;


   //cell position in Angstrom
   double x_cell = program::track_pos_x + cx;
   double y_cell = cy;
   double z_cell = program::track_pos_z  + cz;

   //track sizes/2
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

      //std::cout << x << '\t' << y << '\t' << z << std::endl;

      double Bx = 0.0;
      double By = 0.0;
      double Bz = 0.0;

      for(int k=1; k<=2; k++){

          // predefine power as fixed for loop iteration
          const double m1k = pow(-1,k);

          for(int l=1; l<=2; l++){

             // predefine power as fixed for loop iteration
             const double m1l = pow(-1,l);

             for(int m=1; m<=2; m++){

                const double m1m = pow(-1,m);
                const double m1klm = pow(-1,k+l+m);

                const double xp = x + xb*m1k;
                const double yp = y + yb*m1l;
                const double zp = z + zb*m1m;

                const double xabs = fabs(xp);
                const double yabs = fabs(yp);

                double r = sqrt(xp*xp + yp*yp + zp*zp);

                Bx = Bx + m1klm* log(zp + r);
                By = By + m1klm * sign(yp) * sign(xp) * atan(xabs * zp / (yabs * r));
                Bz = Bz + m1klm* log(xp + r);
                //std::cout <<"BIT\t" << xp << '\t' <<  bit  << '\t' << Bx << '\t' << By << '\t' << Bz << std::endl;


             }
          }
      }
   //std::cout <<"field1:\t" <<  Bx << '\t' <<  By << "\t" << tp::bit_magnetisation[bit] << "\t" << prefactor <<std::endl;
  // std::cout << "\t" << std::endl;
      B[0] = B[0] + Bx*prefactor*tp::bit_magnetisation[bit];
      B[1] = B[1] - By*prefactor*tp::bit_magnetisation[bit];
      B[2] = B[2] + Bz*prefactor*tp::bit_magnetisation[bit];
   //   std::cout <<"field:\t" <<Bx << '\t' << By << '\t' << Bz << '\t' <<  B[0] << '\t' << B[1] << '\t' << B[2] << '\t' << prefactor << '\t' << tp::bit_magnetisation[bit]  << std::endl;

   }

   return B;

}

//------------------------------------------------------------------------------
// Function to set the track field in the local cells at a given step
//------------------------------------------------------------------------------
void set_local_track_field(const int step){
   for (int lc = 0; lc < cells::num_local_cells; lc++){
      const int cell = cells::cell_id_array[lc];
      const double cx = cells::pos_and_mom_array[4*cell+0];
      const double cy = cells::pos_and_mom_array[4*cell+1];
      const double cz = cells::pos_and_mom_array[4*cell+2];
      const std::vector<double> B = calculate_field(cx,cy,cz,step);
      program::track_field_x[cell] = B[0];
      program::track_field_y[cell] = B[1];
      program::track_field_z[cell] = B[2];
   }
}

//------------------------------------------------------------------------------
// Function to write a line of position.txt (root process only): time, track
// position, MR resistance, the track field averaged over all cells with
// magnetic atoms and optionally the field in the last such cell. The cell
// centres are held on every process, so the average does not depend on the
// number of processes.
//------------------------------------------------------------------------------
void write_track_position(std::ofstream& ofile, const int step, const bool last_cell_field){

   if(vmpi::my_rank != 0) return;

   double avB[3] = {0.0, 0.0, 0.0};
   std::vector<double> B(3,0.0);
   int num_cells_with_atoms = 0;

   for (int cell = 0; cell < cells::num_cells; cell++){
      if(cells::num_atoms_in_cell_global[cell] == 0) continue;
      const double cx = cells::pos_and_mom_array[4*cell+0];
      const double cy = cells::pos_and_mom_array[4*cell+1];
      const double cz = cells::pos_and_mom_array[4*cell+2];
      B = calculate_field(cx,cy,cz,step);
      for(int i = 0; i < 3; i++) avB[i] += B[i];
      num_cells_with_atoms++;
   }
   for(int i = 0; i < 3; i++) avB[i] /= double(num_cells_with_atoms);

   ofile << sim::time << '\t' << program::track_pos_z << "\t" << program::track_pos_x << '\t' << micromagnetic::MR_resistance << "\t" << avB[0] << '\t' << avB[1] << '\t' << avB[2];
   if(last_cell_field) ofile << "\t" << B[0] << '\t' << B[1] << '\t' << B[2];
   ofile << std::endl;

}

void tracks(){


     tp::num_bits_per_track = program::track_num_bits_per_track;
     tp::num_tracks = program::track_num_tracks;



     // distance of tracks from read head
     tp::fly_height = program::track_fly_height; // Angstroms

     tp::bit_size =  program::track_bit_size; // size of bits in x-direction (cross track)
    // std::cout << tp::bit_size <<std::endl;
     tp::bit_width = program::track_bit_width; // size of bits in z-direction (down track)
     tp::bit_depth = program::track_bit_depth; // depth of bits along y-direction

     tp::track_gap = program::track_track_gap;
     tp::bit_gap = program::track_bit_gap;

     // centre to centre spacing of bits: across tracks (x) the bit size plus the gap between tracks,
     // down a track (z) the bit width plus the gap between bits
     tp::total_bit_size = tp::bit_size + tp::track_gap;
     tp::total_bit_width = tp::bit_width + tp::bit_gap;
     tp::total_bit_depth = tp::bit_depth;

     // check track parameters (a zero bit size gives a singular field)
     if(tp::num_tracks < 1 || tp::num_bits_per_track < 1){
        err::fatal("Error: disk-tracks needs sim:num-tracks and sim:num-bits-per-track of at least 1. Exiting.");
     }
     if(tp::bit_size <= 0.0 || tp::bit_width <= 0.0 || tp::bit_depth <= 0.0){
        err::fatal("Error: disk-tracks needs sim:track-bit-size-x, -y and -z greater than zero. Exiting.");
     }
     // note: the track field is evaluated at the cell centres in the full size cell arrays, which
     // initialize_modules always builds for this program (dipole::full_cell_fields)

     program::Ms = program::track_Ms;// mu0 Ms in Tesla

     program::track_field_x.resize(cells::num_cells,0.0);
     program::track_field_y.resize(cells::num_cells,0.0);
     program::track_field_z.resize(cells::num_cells,0.0);

   // define total number of bits
    tp::num_bits = (tp::num_bits_per_track)*(tp::num_tracks);


    tp::x_track_pos.resize(tp::num_bits,0.0);
    tp::z_track_pos.resize(tp::num_bits,0.0);
    tp::bit_magnetisation.resize(tp::num_bits,0.0);

    int M = 1;

    std::ifstream ifile;
    ifile.open("track_ms");

    //int i = 0;

    int track_num;
    int bit_num;
    int bits = 0;
    std::vector <double > bitms(tp::num_tracks*tp::num_bits_per_track,0.0);

    if (program::track_ms_file == true){
      zlog << zTs() << "disk-tracks: creating track from file track_ms" << std::endl;

     // read number of atoms in this file
     std::string line; // declare a string to hold line of text
     while(getline(ifile,line) ){
       std::stringstream line_stream(line);
       bits++;
       if (bits > tp::num_bits) {
         err::fatal("Error: track_ms file has more lines than the number of bits (" + std::to_string(tp::num_bits) + "). Exiting.");
       }
       line_stream >> track_num >> bit_num >> M;
   //    std::cout << (track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << '\t' << bitms.size() << "\t" << tp::num_bits_per_track << '\t' <<tp::num_tracks<<std::endl;
       // track and bit numbers in the file count from 1
       if (track_num < 1 || track_num > tp::num_tracks || bit_num < 1 || bit_num > tp::num_bits_per_track){
         err::fatal("Error: track " + std::to_string(track_num) + " bit " + std::to_string(bit_num) +
                    " in track_ms file is out of range (tracks 1-" + std::to_string(tp::num_tracks) +
                    ", bits per track 1-" + std::to_string(tp::num_bits_per_track) + "). Exiting.");
       }
       bitms[(track_num-1)*(tp::num_bits_per_track) + (bit_num-1)] = M;
   //    std::cout <<(track_num-1)*(tp::num_bits_per_track) + (bit_num-1) << std::endl;
     }
   }
   else {
     zlog << zTs() << "disk-tracks: no track file, creating track using alternating Ms values" << std::endl;
     // bits alternate +1, -1 along each track
     for (int track = 0; track < tp::num_tracks; track++){
       for (int bit = 0; bit < tp::num_bits_per_track; bit++){
         bitms[track*tp::num_bits_per_track + bit] = (bit % 2 == 0) ? 1 : -1;
       }
     }
   }

   // bit centres (A) on a grid centred on the origin, tracks along x and bits down each track along z,
   // in the same track-major order as bitms
   for (int track = 0; track < tp::num_tracks; track++){
      for (int b = 0; b < tp::num_bits_per_track; b++){
         const int bit = track*tp::num_bits_per_track + b;
         tp::x_track_pos[bit] = (double(track) - 0.5*double(tp::num_tracks - 1))*tp::total_bit_size;
         tp::z_track_pos[bit] = (double(b) - 0.5*double(tp::num_bits_per_track - 1))*tp::total_bit_width;
         tp::bit_magnetisation[bit] = bitms[bit];
      }
   }
   zlog << zTs() << "disk-tracks: " << tp::num_tracks << " tracks of " << tp::num_bits_per_track << " bits, pitch " <<
           tp::total_bit_size << " A across tracks (x) and " << tp::total_bit_width << " A down track (z), bit centres from (" <<
           tp::x_track_pos[0] << ", " << tp::z_track_pos[0] << ") to (" << tp::x_track_pos[tp::num_bits-1] << ", " <<
           tp::z_track_pos[tp::num_bits-1] << ") A" << std::endl;

    //run a scan over bits
   if (program::LFA == false){
   zlog << zTs() << "disk-tracks: moving track scan over " << cells::num_local_cells << " local cells" << std::endl;

   set_local_track_field(sim::time);

   sim::integrate(sim::equilibration_time);

   std::ofstream ofile;
   if(vmpi::my_rank == 0) ofile.open ("position.txt");
   int step = 0;

   while(sim::time <sim::equilibration_time+sim::total_time){

     program::track_pos_x = program::initial_cross_track_position + program::cross_track_velocity*step;
     program::track_pos_z = program::initial_down_track_position + program::down_track_velocity*step;

    set_local_track_field(step);
    write_track_position(ofile, step, true);

         // Integrate system
  sim::integrate(sim::partial_time);

  //ofile << sim::time << '\t' << program::track_pos_z << "\t" << program::track_pos_x << '\t' << micromagnetic::MR_resistance << "\t" << B[0] << '\t' << B[1] << '\t' << B[2] <<  std::endl;
  step++;

    // Calculate magnetisation statistics
    stats::update();

    // Output data
    vout::data();


   }
 }


//run an LFA simulation
else {
  //   std::cout << "LFAAAAAA" << std::endl;

   std::ofstream ofile;
   if(vmpi::my_rank == 0) ofile.open ("position.txt");

   std::ifstream ifile2;
   ifile2.open("lfa-ms");

   double start_ms = 0;
   double end_ms = 0;
   double ms_step = 0;


   if (ifile2.good()){
     zlog << zTs() << "disk-tracks: creating lfa scan from file" << std::endl;

   sim::integrate(sim::equilibration_time);
    // read number of atoms in this file
    std::string line; // declare a string to hold line of text
    while(getline(ifile2,line) ){
      std::stringstream line_stream(line);

      line_stream >> start_ms >> end_ms >> ms_step;
      //std::cout << "ENTER" << std::endl;
      program::Ms = start_ms;
      ms_step = sqrt(ms_step*ms_step);
      program::LFA_scan_field_step = ms_step;
      //std::cout << start_ms << '\t' << end_ms << "\t" << ms_step <<std::endl;
      if (start_ms > end_ms){
    //  std::cout << "here" <<std::endl;
        while(program::Ms > end_ms){
    //  std::cout << "also here" <<std::endl;

          set_local_track_field(0);
          write_track_position(ofile, 0, false);

          //      std::cout << "and here" <<std::endl;

          // Integrate system
          sim::integrate(sim::partial_time);

          // Calculate magnetisation statistics
          stats::update();

          // Output data
          vout::data();

          program::Ms = program::Ms - program::LFA_scan_field_step;

          }
       }

       else {


         while(program::Ms < end_ms){

          set_local_track_field(0);

          // Integrate system
          sim::integrate(sim::partial_time);

          // Calculate magnetisation statistics
          stats::update();

          // Output data
          vout::data();
          // field written after the step, as before
          write_track_position(ofile, 0, false);

         program::Ms = program::Ms + program::LFA_scan_field_step;

         }
      }

    }
  }
  else {
    err::fatal("Error: disk-tracks with LFA needs an lfa-ms file in the run directory. Exiting.");
  }

  }
  }
}//end of namespace program
