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
      std::string test="laser-FWHM-x";
      std::string test2;
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 5.0e7,"input","0.0 Angstroms - 5 millimetre");
         hamr::internal::fwhm_x = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="laser-FWHM-y";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 5.0e7,"input","0.0 Angstroms - 5 millimetre");
         hamr::internal::fwhm_y = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-speed";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "velocity", 0.0, 1.0e18,"input","0 Angstroms/s - 100000000 m/s");
         hamr::internal::head_speed = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-field-x";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 Angstroms - 1 millimetre");
         hamr::internal::H_bounds_x = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="head-field-y";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 Angstroms - 1 millimetre");
         hamr::internal::H_bounds_y = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="field-rise-time";
      if(word==test){
         double dt = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0,"input","0.01 attosecond - 1 second");
         hamr::internal::H_rise_time = dt;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="field-fall-time";
      if(word==test){
         double dt = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0,"input","0.01 attosecond - 1 second");
         hamr::internal::H_fall_time = dt;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="number-of-bits";
      if(word==test){
         int n = atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(n, word, line, prefix, 0, 1000, "input", "1 - 1000");
         hamr::internal::num_bits = n;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="bit-size";
      test2 = "bit-length";
      if(word==test || word==test2){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::bit_size = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-size";
      test2 = "track-width";
      if(word==test || word==test2){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::track_size = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="track-padding";
      if(word==test){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         hamr::internal::track_padding = f;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="bit-sequence-type";
      if(word==test){
         // single-tone sequence automatically defined given bit_size, track_size and system dimensions
         test="single-tone-predefined";
         if(value==test){
            hamr::internal::create_singletone = true;
            hamr::internal::enabled = true;
            return true;
         }
         test="user-defined";
         if(value==test){
            hamr::internal::create_singletone = false;
            hamr::internal::enabled = true;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - value for \'hamr:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"\"" << std::endl;
            std::cerr << "\t\"single-tone-predefined\"" << std::endl;
            std::cerr << "\t\"user-defined\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="bit-sequence";
      // Accepted values are:
      //  1 -> for bit with polarisation along field direction
      // -1 -> for bit with polarisation opposite to field direction
      //  0 -> for bit where field is not applied
      if(word==test){
         std::vector<int> u;
         // read values from string
         u=vin::integers_from_string(value);
         vin::check_for_valid_bitsequence(u, word, line, prefix, -1, 1, "input", "-1, 0, 1");
         // Store sanitised vector into bit seuqnce
         hamr::internal::bit_sequence.clear();
         hamr::internal::bit_sequence = u;
         hamr::internal::enabled = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="NFT-to-pole-spacing";
      test2 = "NPS";
      if(word==test || word==test2){
         double f = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(f, word, line, prefix, unit, "length", 0.0, 1.0e7,"input","0.0 Angstroms - 1 millimetre");
         hamr::internal::NPS = f;
         hamr::internal::enabled = true;
         return true;
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;
   }
}
