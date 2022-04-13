//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "gpu.hpp"
#include "hamr.hpp"
#include "vio.hpp"

// hamr headers
#include "internal.hpp"

namespace hamr{

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------------------------

   //---------------------------------------------------------------------------------
   // Function for initialising hamr data structures and variables
   //---------------------------------------------------------------------------------
   void initialize(const double Hmin,
                   const double Hmax,
                   const double Tmin,
                   const double Tmax,
                   const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const int num_local_atoms
                  ){

		// If program:hamr-simulation not requested, then stop
      if(!hamr::internal::enabled){
         return;
      }

      // output informative message
      zlog << zTs() << "Initialising data structures for hamr calculation." << std::endl;

      // check for prior initialisation
      if(hamr::internal::initialised){
         zlog << zTs() << "Warning: hamr module already initialised. Continuing." << std::endl;
         return;
      }

      //-------------------------------------------------------
      // Save value of local num atoms and resize field arrays
      //-------------------------------------------------------
      hamr::internal::num_local_atoms = num_local_atoms;
      hamr::internal::x_field_array.resize(num_local_atoms, 0.0); // arrays to store atomic hamr fields
      hamr::internal::y_field_array.resize(num_local_atoms, 0.0);
      hamr::internal::z_field_array.resize(num_local_atoms, 0.0);

      hamr::internal::system_dimensions_x = system_dimensions_x;
      hamr::internal::system_dimensions_y = system_dimensions_y;
      hamr::internal::system_dimensions_z = system_dimensions_z;
      hamr::internal::atom_coords_x = atom_coords_x;
      hamr::internal::atom_coords_y = atom_coords_y;
      hamr::internal::atom_coords_z = atom_coords_z;

      hamr::internal::atom_type_array = atom_type_array;
      hamr::internal::Hmin = Hmin;
      hamr::internal::Hmax = Hmax;
      hamr::internal::Tmin = Tmin;
      hamr::internal::Tmax = Tmax;

      // Calibrate head region to be not larger than system size
      if(hamr::internal::H_bounds_x > hamr::internal::system_dimensions_x){ hamr::internal::H_bounds_x = hamr::internal::system_dimensions_x;}
      if(hamr::internal::H_bounds_y > hamr::internal::system_dimensions_y){ hamr::internal::H_bounds_y = hamr::internal::system_dimensions_y;}

      const double one_over_sqrt = 1.0/sqrt(8.0*log(2.0));
      hamr::internal::laser_sigma_x = hamr::internal::fwhm_x * one_over_sqrt;
      hamr::internal::laser_sigma_y = hamr::internal::fwhm_y * one_over_sqrt;

      // Calculate max number of allowed tracks and bits-per-track in the system
      hamr::internal::num_tracks = int(floor((hamr::internal::system_dimensions_y - 1.0 - 2.0*hamr::internal::track_padding)/(hamr::internal::track_size)));
      hamr::internal::bits_per_track = int(floor((hamr::internal::system_dimensions_x - 1.0)/hamr::internal::bit_size));

      // If number of tracks determined by the requested number of bits is less than the maximum number of tracks allowed,
      // update the number of tracks to be written
      const int num_tracks_try = int(round(1.0*hamr::internal::num_bits/(1.0*hamr::internal::bits_per_track)));
      if(num_tracks_try < hamr::internal::num_tracks){
         std::cout << "Warning: Reducing number of track(s) from " << hamr::internal::num_tracks << " to " << num_tracks_try << std::endl;
         zlog << zTs() << "Warning: Reducing number of track(s) from " << hamr::internal::num_tracks << " to " << num_tracks_try << std::endl;
         hamr::internal::num_tracks = num_tracks_try;
      }

      // Check if need to create singletone
      hamr::internal::create_singletone_vector();

      // Check that provided bit sequence is consistent with input data and system dimensions
      hamr::internal::check_sequence_length();

      // Set initialised flag
      hamr::internal::initialised = true;

      return;
   }

} // end of hamr namespace
