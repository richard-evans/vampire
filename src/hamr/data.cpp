//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "hamr.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------------------------
   bool head_laser_on   = false;

   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for hamr calculation
      //-----------------------------------------------------------------------------
      bool enabled = false;
      bool initialised = false;
      bool create_singletone = false;

      int num_bits = 0;
      int bits_per_track = 0;
      int num_tracks = 0;

      double bit_size = 0.0;
      double track_size = 0.0;
      double head_position_x = 0.0;
      double head_position_y = 0.0;
      double head_speed = 30.0e10; // A/s
      double fwhm_x = 200.0; // A
      double fwhm_y = 200.0; // A
      double laser_sigma_x = fwhm_x/sqrt(8.0*log(2.0));
      double laser_sigma_y = fwhm_y/sqrt(8.0*log(2.0));
      double H_bounds_x = 200.0; //
      double H_bounds_y = 200.0; //
      double H_rise_time = 1.0e-12;
      double H_fall_time = 1.0e-12;
      double NPS = 0.0; // NFT to pole spacing
      double track_padding = 0.0;
      std::vector<int> bit_sequence;

      int num_local_atoms;
      double Tmin = 0.0;
      double Tmax = 0.0;
      double Hmin = 0.0;
      double Hmax = 0.0;
      double system_dimensions_x;
      double system_dimensions_y;
      double system_dimensions_z;
      std::vector<double> x_field_array; // arrays to store atomic spin torque field
      std::vector<double> y_field_array;
      std::vector<double> z_field_array;
      std::vector<int> atom_type_array;
      std::vector<double> atom_coords_x;
      std::vector<double> atom_coords_y;
      std::vector<double> atom_coords_z;

   } // end of internal namespace
} // end of hamr namespace
