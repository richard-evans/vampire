//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "program.hpp"
#include "errors.hpp"
#include "vio.hpp"

// program module headers
#include "internal.hpp"

namespace program{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for program module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false - note use of sim: module name for semantic reasons
      // program variables are generally related to simulations but a logically separated for ease of the programmer
      std::string prefix="sim";
      if(key!=prefix) return false;
      //-------------------------------------------------------------------
      std::string test="program";
      if(word==test){
         test="benchmark";
         if(value==test){
            program::program=0;
            return true;
         }
         test="time-series";
         if(value==test){
            program::program=1;
            return true;
         }
         test="hysteresis-loop";
         if(value==test){
            program::program=2;
            return true;
         }
         test="static-hysteresis-loop";
         if(value==test){
            program::program=3;
            return true;
         }
         test="curie-temperature";
         if(value==test){
            program::program=4;
            return true;
         }
         test="field-cool";
         if(value==test){
            program::program=5;
            return true;
         }
         test="localised-field-cool";
         if(value==test){
            program::program=16;
            return true;
         }
         test="laser-pulse";
         if(value==test){
            program::program=6;
            return true;
         }
         test="hamr-simulation";
         if(value==test){
            program::program=7;
            return true;
         }
         test="cmc-anisotropy";
         if(value==test){
            program::program=8;
            return true;
         }
         test="hybrid-cmc";
         if(value==test){
            program::program=9;
            return true;
         }
         test="reverse-hybrid-cmc";
         if(value==test){
            program::program=10;
            return true;
         }
         test="LaGrange-Multiplier";
         if(value==test){
            program::program=11;
            return true;
         }
         test="partial-hysteresis-loop";
         if(value==test){
            program::program=12;
            return true;
         }
         test="localised-temperature-pulse";
         if(value==test){
            program::program=13;
            return true;
         }
         test="effective-damping";
         if(value==test){
            program::program=14;
            return true;
         }
         test="fmr";
         if(value==test){
            program::program=15;
            return true;
         }
         test="diagnostic-boltzmann";
         if(value==test){
            program::program=50;
            return true;
         }
         test="setting";
         if(value==test){
            program::program=51;
            return true;
         }
         test="domain-wall";
         if(value==test){
            program::program=52;
            return true;
         }
         test="exchange-stiffness";
         if(value==test){
            program::program=53;
            return true;
         }
         test="electrical-pulse";
         if(value==test){
            program::program = 17;
            return true;
         }
         test="mm-A-calculation";
         if(value==test){
            program::program=54;
            return true;
         }
         test="field-sweep";
         if(value==test){
            program::program=70;
            return true;
         }
         test="disk-tracks";
         if(value==test){
            program::program=72;
            return true;
         }
         test="diagnostic-boltzmann-micromagnetic-llg";
         if(value==test){
            program::program=73;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cout << word << '\t' << test << std::endl;
            std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"benchmark\"" << std::endl;
            std::cerr << "\t\"cmc-anisotropy\"" << std::endl;
            std::cerr << "\t\"curie-temperature\"" << std::endl;
            std::cerr << "\t\"domain-wall\"" << std::endl;
            std::cerr << "\t\"effective-damping\"" << std::endl;
            std::cerr << "\t\"electrical-pulse\"" << std::endl;
            std::cerr << "\t\"exchange-stiffness\"" << std::endl;
            std::cerr << "\t\"field-cool\"" << std::endl;
            std::cerr << "\t\"laser-pulse\"" << std::endl;
            std::cerr << "\t\"localised-field-cool\"" << std::endl;
            std::cerr << "\t\"localised-temperature-pulse\"" << std::endl;
            std::cerr << "\t\"time-series\"" << std::endl;
            std::cerr << "\t\"hysteresis-loop\"" << std::endl;
            std::cerr << "\t\"partial-hysteresis-loop\"" << std::endl;
            std::cerr << "\t\"hybrid-cmc\"" << std::endl;
            std::cerr << "\t\"reverse-hybrid-cmc\"" << std::endl;
            std::cerr << "\t\"static-hysteresis-loop\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-time";
      if(word == test){
         double pt = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(pt, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_time = pt;
         return true;
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-rise-time";
      if(word == test){
         double rt = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(rt, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_rise_time = rt;
         return true;
      }
      //-------------------------------------------------------------------
      test = "electrical-pulse-fall-time";
      if(word == test){
         double ft = atof(value.c_str()); // convert string to uint64_t
         // Test for valid range
         vin::check_for_valid_positive_value(ft, word, line, prefix, unit, "time", 0.0, 1.0,"input","0.0 - 1s");
         // save sanitized value
         program::internal::electrical_pulse_fall_time = ft;
         return true;
      }
      //--------------------------------------------------------------------
      test="exchange-stiffness-maximum-angle";
      if(word==test){
         double ma = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(ma, word, line, prefix, unit, "", 0.0, 180.1,"input","0 - 180 degrees");
         program::internal::exchange_stiffness_max_constraint_angle = ma;
         return true;
      }
      //--------------------------------------------------------------------
      test="exchange-stiffness-angle-increment";
      if(word==test){
         double ai = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(ai, word, line, prefix, unit, "", 1.0, 90.0,"input","1 - 90 degrees");
         program::internal::exchange_stiffness_delta_constraint_angle = ai;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of program namespace
