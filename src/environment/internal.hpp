//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ENVIRONMENT_INTERNAL_H_
#define ENVIRONMENT_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the environment module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <fstream>
#include <iostream>
#ifdef FFT
#include <fftw3.h>
#endif

// Vampire headers
#include "environment.hpp"

// environment module headers
#include "internal.hpp"

using namespace std;

namespace environment{

   namespace internal{

     extern int num_shields;

     extern std::vector <string> shield_shape;
     extern std::vector <double> shield_ms;
     extern std::vector <double> shield_Tc;
     extern std::vector< std::vector <double> > shield_A;
     extern std::vector <double> shield_ku;
     extern std::vector <double> shield_alpha;
     extern std::vector <double> shield_gamma;
     extern std::vector <double> shield_max_x;
     extern std::vector <double> shield_max_y;
     extern std::vector <double> shield_max_z;
     extern std::vector <double> shield_min_x;
     extern std::vector <double> shield_min_y;
     extern std::vector <double> shield_min_z;
     extern std::vector <double> shield_Hext_x;
     extern std::vector <double> shield_Hext_y;
     extern std::vector <double> shield_Hext_z;
     extern std::vector <double> shield_max_cell_size;
     extern std::vector <double> shield_min_cell_size;
     extern std::vector <string> pos_or_neg;
     extern std::vector <int> shield_number;
     extern std::vector <int> H_strength;



      extern std::vector < int > cell_dx;
      extern std::vector < int > cell_dy;
      extern std::vector < int > cell_dz;

      extern std::vector < std::vector < std::vector<int> > > idarray;

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      extern std::vector < double > dim;
      extern std::vector < double > cell_size_x;
      extern std::vector < double > cell_size_y;
      extern std::vector < double > cell_size_z;
      extern std::vector < double > initial_spin_x;
      extern std::vector < double > initial_spin_y;
      extern std::vector < double > initial_spin_z;
      extern std::vector < double > ext_field;


      extern std::vector < double > dipole_field_x;
      extern std::vector < double > dipole_field_y;
      extern std::vector < double > dipole_field_z;

      extern std::vector < double > bias_field_x;
      extern std::vector < double > bias_field_y;
      extern std::vector < double > bias_field_z;
      //3d vector to store the shift in atomistic position from 0,0,0 with respect to the environment module
      extern std::vector < double > shift;

      extern bool square_shields;
      extern bool expoential_shields;
      extern int gap;

      extern bool LFA_scan;
      extern double eightPI_three_cell_volume;
      extern std::vector < double > cell_volume;
      extern double env_field;
      extern std::vector < double > env_field_uv;
      //if initial spins are random
      extern std::vector < bool >  random_spins;

      //stores cell properties
      extern int num_cells;
      extern double exchange_constant;
      extern std::vector < double > ku;// -1e-23;
      extern std::vector < double > Ms;// = 1e-21;

      extern double m_e;
      extern double gamma;
      extern double alpha;
      extern int num_env_cells;
      extern std::vector < double > one_o_chi_para;
      extern std::vector < double > one_o_chi_perp;
      extern double alpha_para;
      extern double alpha_perp;
      extern double cell_size;
      extern bool env_output_info;
      extern std::vector <int> list_of_mm_cells_in_env;
      //array to store mag
      extern std::vector < double > x_mag_array;
      extern std::vector < double > y_mag_array;
      extern std::vector < double > z_mag_array;
      #ifdef FFT
      //FT for magnetisation
      extern fftw_plan MxP,MyP,MzP;
      // performs the backward transform to give the dipole field, Hx, Hy, Hz
      extern fftw_plan HxP,HyP,HzP;
      #endif
//      extern fft_cell_id_array;

      extern std::vector<double> cell_coords_array_x; /// arrays to store cells positions
      extern std::vector<double> cell_coords_array_y;
      extern std::vector<double> cell_coords_array_z;

      //stores nearest neighbours
      extern std::vector<double> neighbour_list_start_index;
      extern std::vector<double> neighbour_list_end_index;
      extern std::vector<double> neighbour_list_array;

      extern ofstream o_file;
      #ifdef FFT
      extern fftw_complex *Mx_in; //3D Array for m_in
      extern fftw_complex *My_in;
      extern fftw_complex *Mz_in;

      extern fftw_complex *Hx_in; //3D Array for Hin
      extern fftw_complex *Hy_in;
      extern fftw_complex *Hz_in;

      extern fftw_complex *Mx_out; //3D Array for mout
      extern fftw_complex *My_out;
      extern fftw_complex *Mz_out;

      extern fftw_complex *Hx_out; //3D Array for Hout
      extern fftw_complex *Hy_out;
      extern fftw_complex *Hz_out;


      extern fftw_complex *N2xx0; //3D Array for dipolar field
      extern fftw_complex *N2xy0;
      extern fftw_complex *N2xz0;

      extern fftw_complex *N2yx0; //3D Array for dipolar field
      extern fftw_complex *N2yy0;
      extern fftw_complex *N2yz0;

      extern fftw_complex *N2zx0; //3D Array for dipolar field
      extern fftw_complex *N2zy0;
      extern fftw_complex *N2zz0;

      extern fftw_complex *N2xx; //3D Array for dipolar field
      extern fftw_complex *N2xy;
      extern fftw_complex *N2xz;

      extern fftw_complex *N2yx; //3D Array for dipolar field
      extern fftw_complex *N2yy;
      extern fftw_complex *N2yz;

      extern fftw_complex *N2zx; //3D Array for dipolar field
      extern fftw_complex *N2zy;
      extern fftw_complex *N2zz;

      #endif
      extern int num_cells_x;
      extern int num_cells_y;
      extern int num_cells_z;

      extern int eight_num_cells;

      extern std::vector < double > rij_tensor_xx;
      extern std::vector < double > rij_tensor_xy;
      extern std::vector < double > rij_tensor_xz;

      extern std::vector < double > rij_tensor_yy;
      extern std::vector < double > rij_tensor_yz;
      extern std::vector < double > rij_tensor_zz;


      extern std::vector < int > none_atomistic_cells;
      extern std::vector < int > atomistic_cells;

      extern std::vector < int > list_env_cell_atomistic_cell;
      extern std::vector < int > env_cell_is_in_atomistic_region;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      double calculate_chi_perp(double T, int cell);
      double calculate_chi_para(double T, int cell);
      int initialise_demag_fields();
      int calculate_demag_fields();
      int output();
      int read_in_shield_info();

      int in_shield(double x, double y, double z, int shield);
      int bias_shields();

      std::vector<double> calculate_llb_fields(std::vector <double>& m,
                                               double temperature,
                                               int cell,
                                               std::vector<double>& x_array,
                                               std::vector<double>& y_array,
                                               std::vector<double>& z_array);

      std::vector<double> calculate_field_env(int celli, int cellj);
      std::vector<double> calculate_field_mm(int celli, int cellj);
      std::vector < std::vector < double> > calculate_corners(double x, double y, double z, double cell_size_x, double cell_size_y, double cell_size_z);

      } // end of internal namespace

   } // end of environment namespace

   #endif //ENVIRONMENT_INTERNAL_H_
