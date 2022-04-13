//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "hamr.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace hamr{

   //--------------------------------------------------------------------------------
   // Function to get head position along x-direction (downtrack)
   //--------------------------------------------------------------------------------
   double get_head_position_x(){
      return internal::head_position_x;
   }

   //--------------------------------------------------------------------------------
   // Function to get head position along y-direction (crosstrack)
   //--------------------------------------------------------------------------------
   double get_head_position_y(){
      return internal::head_position_y;
   }

   //--------------------------------------------------------------------------------
   // Function to get head speed
   //--------------------------------------------------------------------------------
   double get_head_speed(){
      return internal::head_speed;
   }

   //--------------------------------------------------------------------------------
   // Function to get FWHM in x-direction
   //--------------------------------------------------------------------------------
   double get_FWHM_x(){
      return internal::fwhm_x;
   }

   //--------------------------------------------------------------------------------
   // Function to get FWHM in y-direction
   //--------------------------------------------------------------------------------
   double get_FWHM_y(){
      return internal::fwhm_y;
   }

   //--------------------------------------------------------------------------------
   // Function to get standard devieation of temperature profile in x-direction
   //--------------------------------------------------------------------------------
   double get_laser_sigma_x(){
      return internal::laser_sigma_x;
   }

   //--------------------------------------------------------------------------------
   // Function to get standard devieation of temperature profile in y-direction
   //--------------------------------------------------------------------------------
   double get_laser_sigma_y(){
      return internal::laser_sigma_y;
   }

   //--------------------------------------------------------------------------------
   // Function to get region of application of external field in x-direction
   //--------------------------------------------------------------------------------
   double get_field_bounds_x(){
      return internal::H_bounds_x;
   }

   //--------------------------------------------------------------------------------
   // Function to get region of application of external field in x-direction
   //--------------------------------------------------------------------------------
   double get_field_bounds_y(){
      return internal::H_bounds_y;
   }

   //--------------------------------------------------------------------------------
   // Function to get bit length (downtrack)
   //--------------------------------------------------------------------------------
   double get_bit_size(){
      return internal::bit_size;
   }

   //--------------------------------------------------------------------------------
   // Function to get track width (crosstrack)
   //--------------------------------------------------------------------------------
   double get_track_size(){
      return internal::track_size;
   }

   //--------------------------------------------------------------------------------
   // Function to get external field ramp time
   //--------------------------------------------------------------------------------
   double get_field_rise_time(){
      return internal::H_rise_time;
   }

   //--------------------------------------------------------------------------------
   // Function to get external field ramp time
   //--------------------------------------------------------------------------------
   double get_field_fall_time(){
      return internal::H_fall_time;
   }

   //--------------------------------------------------------------------------------
   // Function to get external field ramp time
   //--------------------------------------------------------------------------------
   double get_NPS(){
      return internal::NPS;
   }

   //--------------------------------------------------------------------------------
   // Function to get track padding
   //--------------------------------------------------------------------------------
   double get_track_padding(){
      return internal::track_padding;
   }

   //--------------------------------------------------------------------------------
   // Function to get hamr initialised state on cpu
   //--------------------------------------------------------------------------------
   bool get_initialisation_state(){
      return hamr::internal::initialised;
   }

} // end of hamr namespace
