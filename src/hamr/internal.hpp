#ifndef HAMR_INTERNAL_H_
#define HAMR_INTERNAL_H_
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

   namespace internal{
      //-----------------------------------------------------------------------------
      // Function to calculate Gaussian temperature profile
      //-----------------------------------------------------------------------------
      double calculate_gaussian_profile(const int atom, 
                                       const double Tmin, 
                                       const double DeltaT); 

      //-----------------------------------------------------------------------------
      // Function to calculate the external field with trapezoidal temporal profile
      //-----------------------------------------------------------------------------
      void update_field_time_trapz_profile(const uint64_t field_time,   
                                          const uint64_t ramp_time,
                                          const uint64_t ramp_time_end,
                                          const uint64_t pre_write_time,
                                          const uint64_t write_time,
                                          double &H_applied);

      //-----------------------------------------------------------------------------
      // Function to apply field only to atoms within box around centre of writer coil
      //-----------------------------------------------------------------------------
      void apply_field_spatial_box(const int start_index,
					                  const int end_index,
											const double Hloc_parity_field,
											const double Hvecx,
											const double Hvecy,
											const double Hvecz,
					                  std::vector<double>& x_total_external_field_array,
					                  std::vector<double>& y_total_external_field_array,
					                  std::vector<double>& z_total_external_field_array);

      //-----------------------------------------------------------------------------
      // Shared variables used for hamr calculation
      //-----------------------------------------------------------------------------
      extern bool initialised; // flag set if initialised
      extern bool single_bit;

      extern double head_position[2];
      extern double head_speed;
      extern double laser_peak_time;
      extern double fwhm_x;
      extern double fwhm_y;
      extern double H_bounds_x;
      extern double H_bounds_y;
      extern double H_osc_amplit;
      extern double H_ramp_time;
      extern double bit_spacing_x;
      extern double bit_spacing_y;
      extern double Hmin;
      extern double Hmax;

      extern int num_local_atoms;

      extern std::vector<double> x_field_array; /// arrays to store atomic hamr field
      extern std::vector<double> y_field_array;
      extern std::vector<double> z_field_array;

      extern double system_dimensions[3];
      // extern std::vector<double> system_dimensions_x; /// 
      // extern std::vector<double> system_dimensions_y;
      // extern std::vector<double> system_dimensions_z;

      extern std::vector<double> atom_coords_x; /// 
      extern std::vector<double> atom_coords_y;
      extern std::vector<double> atom_coords_z;

      extern std::vector<int> atom_type_array;


   } // end of internal namespace
} // end of hamr namespace

#endif //HAMR_INTERNAL_H_