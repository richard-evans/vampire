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

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "environment.hpp"

// environment module headers
#include "internal.hpp"

namespace environment{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   std::vector < double > environment_field_x;
   std::vector < double > environment_field_y;
   std::vector < double > environment_field_z;

   std::vector < double > atomic_field_x;
   std::vector < double > atomic_field_y;
   std::vector < double > atomic_field_z;

   int demag_update_rate = 1000;
   bool enabled = false;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside environment module
      //------------------------------------------------------------------------

      std::vector < double > dim(3,6.0);
      std::vector < double > cell_size(3,2.0);

      std::vector < double > initial_spin(3,0.0);

      std::vector < double > ext_field(3,0.0);

      std::vector < double > dipole_field_x;
      std::vector < double > dipole_field_y;
      std::vector < double > dipole_field_z;

      std::vector < double > shift(3,0.0);
      double eightPI_three_cell_volume;
      double cell_volume;

      bool random_spins = true;

      double one_o_chi_para;
      double one_o_chi_perp;

      int num_cells;
      double A = -180;
      double ku= -1e-23;
      double Tc = 600;
      double Ms = 1e-21;
      double gamma = 1.0;
      double alpha = 1.0;
      int num_env_cells = 0.0;
      double m_e;
      double exchange_constant;
      double alpha_para;
      double alpha_perp;

      std::vector < double > x_mag_array;
      std::vector < double > y_mag_array;
      std::vector < double > z_mag_array;

      std::vector<double> cell_coords_array_x; /// arrays to store cells positions
      std::vector<double> cell_coords_array_y;
      std::vector<double> cell_coords_array_z;

      std::vector<double> neighbour_list_start_index;
      std::vector<double> neighbour_list_end_index;
      std::vector<double> neighbour_list_array;

      ofstream o_file;


      Array3D<fftw_complex> Nxx0; //3D Array for dipolar field
      Array3D<fftw_complex> Nxy0;
      Array3D<fftw_complex> Nxz0;

      Array3D<fftw_complex> Nyx0; //3D Array for dipolar field
      Array3D<fftw_complex> Nyy0;
      Array3D<fftw_complex> Nyz0;

      Array3D<fftw_complex> Nzx0; //3D Array for dipolar field
      Array3D<fftw_complex> Nzy0;
      Array3D<fftw_complex> Nzz0;

      Array3D<fftw_complex> Nxx; //3D Array for dipolar field
      Array3D<fftw_complex> Nxy;
      Array3D<fftw_complex> Nxz;

      Array3D<fftw_complex> Nyx; //3D Array for dipolar field
      Array3D<fftw_complex> Nyy;
      Array3D<fftw_complex> Nyz;

      Array3D<fftw_complex> Nzx; //3D Array for dipolar field
      Array3D<fftw_complex> Nzy;
      Array3D<fftw_complex> Nzz;

      int num_cells_x;
      int num_cells_y;
      int num_cells_z;

      int eight_num_cells;


      std::vector < int > none_atomistic_cells;
      std::vector < int > atomistic_cells;
      std::vector < int > list_env_cell_atomistic_cell;
      std::vector < int > env_cell_is_in_atomistic_region;


   } // end of internal namespace

} // end of environment namespace
