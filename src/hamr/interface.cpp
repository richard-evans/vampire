//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "hamr.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace hamr{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="hamr";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      //--------------------------------------------------------------------
      std::string test="laser-peak-time";
      if(word==test){
         double dt = atof(value.c_str());
         // Test for valid range                                                                                                                            
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
         hamr::internal::laser_peak_time = dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="laser-FWHM-x";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 5.0e7,"input","0.1 Angstroms - 5 millimetre");
         hamr::internal::fwhm_x = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="laser-FWHM-y";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 5.0e7,"input","0.1 Angstroms - 5 millimetre");
         hamr::internal::fwhm_y = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-speed";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "velocity", 0.0, 1.0e18,"input","0 Angstroms/s - 100000000 m/s");
         hamr::internal::head_speed = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-field-x";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::H_bounds_x = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-field-y";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::H_bounds_y = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="field-oscillation-frequency";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::H_osc_amplit = f;
         return true;
      }
      //--------------------------------------------------------------------
      test="field-ramp-time";
      if(word==test){
         double dt = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
         hamr::internal::H_ramp_time = dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="bit-spacing-x";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::bit_spacing_x = f;
         return true;                                                                                                                                       
      }
      //--------------------------------------------------------------------
      test="bit-spacing-y";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::bit_spacing_y = f;
         return true;                                                                                                                                       
      }
      //--------------------------------------------------------------------
      test="single-bit";
      if(word==test){
         // hamr::internal::single_bit = true;
         // return true;
         bool tf = vin::check_for_valid_bool(value, word, line, prefix, "hamr");
         hamr::internal::single_bit = tf;
         return true;
      }
      // //--------------------------------------------------------------------
      // test="continuous";
      // if(word==test){
      //    // hamr::continuous = true;
      //    // return true;
      //    bool tf = vin::check_for_valid_bool(value, word, line, prefix, "hamr");
      //    hamr::continuous = tf;
      //    return true;
      // }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;
   }
}