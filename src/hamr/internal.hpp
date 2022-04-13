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

      void apply_temperature_profile(const int start_index, const int end_index, const double Tmin, const double DeltaT);

      //-----------------------------------------------------------------------------
      // Function to calculate the external field with trapezoidal temporal profile
      //-----------------------------------------------------------------------------
      double update_field_time_trapz_profile(const uint64_t current_time, const uint64_t rise_time, const uint64_t fall_time, const uint64_t bit_time);

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
      // Function to create singletone bit sequence
      //-----------------------------------------------------------------------------
		void create_singletone_vector();

      //-----------------------------------------------------------------------------
      // Function to check that user defined bit sequence is consistent with system
		// size and "number-of-bits" parameter
      //-----------------------------------------------------------------------------
		void check_sequence_length();

      //-----------------------------------------------------------------------------
      // Shared variables used for hamr calculation
      //-----------------------------------------------------------------------------
      extern bool enabled; /// flag set if initialised
      extern bool initialised; /// flag set if initialised
      extern bool create_singletone;  /// flag to create single tone

      extern int num_bits; /// total number of bits
      extern int bits_per_track; /// number of bits per track
      extern int num_tracks; /// number of tracks

      extern double bit_size; /// bit size
      extern double track_size; /// track size
      extern double head_position_x; /// head centre coordinate in x
      extern double head_position_y; /// head centre coordinate in x
      extern double head_speed; /// linear sweeping speed in Ang/s
      extern double fwhm_x;  /// fwhm of Gaussian distribution
      extern double fwhm_y;  /// fwhm of Gaussian distribution
      extern double laser_sigma_x;  /// standard deviation of Gaussian distribution from fwhm/sqrt(8*log(2))
      extern double laser_sigma_y;  /// standard deviation of Gaussian distribution from fwhm/sqrt(8*log(2))
      extern double H_bounds_x;  /// Dimension in x-direction of box where magnetic field is applied
      extern double H_bounds_y;  /// Dimension in y-direction of box where magnetic field is applied
      extern double H_rise_time;  /// Rise time in trapezoidal time profile of external field
      extern double H_fall_time;  /// Fall time in trapezoidal time profile of external field
      extern double NPS;  /// NFT to pole spacing
      extern double track_padding;  ///spacing between edges of the system in y-direction and track(s)
      extern std::vector<int> bit_sequence;  /// vector to store bits

      extern int num_local_atoms;
      extern double Tmin;   /// Minimum temperature of thermal profile
      extern double Tmax;   /// Maximum temperature of thermal profile
      extern double Hmin;   /// Minimum value of field
      extern double Hmax;   /// Maximum value of field
      extern double system_dimensions_x; ///
      extern double system_dimensions_y;
      extern double system_dimensions_z;
      extern std::vector<double> atom_coords_x; ///
      extern std::vector<double> atom_coords_y;
      extern std::vector<double> atom_coords_z;
      extern std::vector<int> atom_type_array;
      extern std::vector<double> x_field_array; /// arrays to store atomic hamr field
      extern std::vector<double> y_field_array;
      extern std::vector<double> z_field_array;

   } // end of internal namespace
} // end of hamr namespace

#endif //HAMR_INTERNAL_H_
