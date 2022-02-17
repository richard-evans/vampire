//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

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
      bool initialised = false;
      bool single_bit = false;
      double head_position[2] = {0.0, system_dimensions[1]*0.5}; // A
      double head_speed = 30.0; // m/s
      double laser_peak_time = 500.0e-12; //s
      double fwhm_x = 200.0; // A
      double fwhm_y = 200.0; // A
      double H_bounds_x = system_dimensions[0]; // 
      double H_bounds_y = system_dimensions[1]; // 
      double H_osc_amplit = H_bounds_x; // A
      double H_ramp_time = 1.0e-12;
      double bit_spacing_x = 0.0;
      double bit_spacing_y = 0.0;

      int num_local_atoms;
      double system_dimensions[3];
      double Hmin = 0.0;
      double Hmax = 0.0;

      std::vector<double> x_field_array; // arrays to store atomic spin torque field
      std::vector<double> y_field_array;
      std::vector<double> z_field_array;

      std::vector<int> atom_type_array;
      std::vector<double> atom_coords_x;
      std::vector<double> atom_coords_y;
      std::vector<double> atom_coords_z;

   } // end of internal namespace
} // end of hamr namespace

