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
#include "array3d.h"


// Vampire headers
#include "environment.hpp"

// environment module headers
#include "internal.hpp"

using namespace std;

namespace environment{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      extern std::vector < double > dim;
      extern std::vector < double > cell_size;
      extern std::vector < double > initial_spin;
      extern std::vector < double > ext_field;


      extern std::vector < double > dipole_field_x;
      extern std::vector < double > dipole_field_y;
      extern std::vector < double > dipole_field_z;

      extern std::vector < double > shift;

      extern double eightPI_three_cell_volume;
      extern double cell_volume;

      extern bool random_spins;

      extern int num_cells;
      extern double A;
      extern double exchange_constant;
      extern double ku;
      extern double Tc;
      extern double Ms;

      extern double m_e;

      extern double gamma;
      extern double alpha;

      extern int num_env_cells;

      extern double one_o_chi_para;
      extern double one_o_chi_perp;

      extern double alpha_para;
      extern double alpha_perp;

      extern std::vector < double > x_mag_array;
      extern std::vector < double > y_mag_array;
      extern std::vector < double > z_mag_array;


      extern std::vector<double> cell_coords_array_x; /// arrays to store cells positions
      extern std::vector<double> cell_coords_array_y;
      extern std::vector<double> cell_coords_array_z;

      extern std::vector<double> neighbour_list_start_index;
      extern std::vector<double> neighbour_list_end_index;
      extern std::vector<double> neighbour_list_array;

      extern ofstream o_file;


      extern Array3D<fftw_complex> Nxx0; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nxy0;
      extern Array3D<fftw_complex> Nxz0;

      extern Array3D<fftw_complex> Nyx0; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nyy0;
      extern Array3D<fftw_complex> Nyz0;

      extern Array3D<fftw_complex> Nzx0; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nzy0;
      extern Array3D<fftw_complex> Nzz0;

      extern Array3D<fftw_complex> Nxx; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nxy;
      extern Array3D<fftw_complex> Nxz;

      extern Array3D<fftw_complex> Nyx; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nyy;
      extern Array3D<fftw_complex> Nyz;

      extern Array3D<fftw_complex> Nzx; //3D Array for dipolar field
      extern Array3D<fftw_complex> Nzy;
      extern Array3D<fftw_complex> Nzz;

      extern int num_cells_x;
      extern int num_cells_y;
      extern int num_cells_z;

      extern int eight_num_cells;


      extern std::vector < int > none_atomistic_cells;
      extern std::vector < int > atomistic_cells;

      extern std::vector < int > list_env_cell_atomistic_cell;
      extern std::vector < int > env_cell_is_in_atomistic_region;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      double calculate_chi_perp(double T);
      double calculate_chi_para(double T);
      int initialise_demag_fields();
      int calculate_demag_fields();
      int output();


      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                                double temperature,
                                                int cell,
                                                std::vector<double> x_array,
                                                std::vector<double> y_array,
                                                std::vector<double> z_array);

   } // end of internal namespace

} // end of environment namespace

#endif //ENVIRONMENT_INTERNAL_H_
