//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <string>

// Vampire headers
#include "errors.hpp"
#include "gpu.hpp"
#include "vio.hpp"

namespace gpu{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ltmp settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="gpu";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="acceleration";
      if(word==test){
         test="on";
         // enable (default for gpu compiled version)
         if(value==test){
            gpu::acceleration = true;
            return true;
         }
         test="off";
         // disable (default for cpu compiled version)
         if(value==test){
            gpu::acceleration = false;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"on\"" << std::endl;
            std::cerr << "\t\"off\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      //--------------------------------------------------------------------
      test="calculate-statistics-on-cpu";
      if(word==test){
         test="";
         if(value==test){
            gpu::cpu_stats = true;
            return true;
         }
         test="true";
         if(value==test){
            gpu::cpu_stats = true;
            return true;
         }
         // default
         test="false";
         if(value==test){
            gpu::cpu_stats = false;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"true\"" << std::endl;
            std::cerr << "\t\"false\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      test = "num-threads";
      if (word == test) {
         int val = atoi(value.c_str());
         if (val > 0) {
            gpu::num_threads = val;
            return true;
         } else {
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ':' << word << "\' must be a positive integer" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      test = "num-threads-other";
      if (word == test) {
         int val = atoi(value.c_str());
         if (val > 0) {
            gpu::num_threads_other = val;
            return true;
         } else {
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ':' << word << "\' must be a positive integer" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      test = "platform";
      if (word == test) {
         int val = atoi(value.c_str());
         if (val >= 0) {
            gpu::platform = val;
            return true;
         } else {
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ':' << word << "\' must be zero or a positive integer" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      //--------------------------------------------------------------------
      test="device";
      if(word==test){
         int dev=atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(dev, word, line, prefix, 0, 16,"input","0-16");
         gpu::device = dev;
         return true;
      }

      test="platform-other";
      if (word == test) {
         int val = atoi(value.c_str());
         if (val >= 0) {
            gpu::platform_other = val;
            return true;
         } else {
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ':' << word << "\' must be zero or a positive integer" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }

      test="device-other";
      if(word==test){
         int dev=atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(dev, word, line, prefix, 0, 16,"input","0-16");
         gpu::device_other = dev;
         return true;
      }

      //--------------------------------------------------------------------
      else{
         terminaltextcolor(RED);
         std::cerr << "Error - Unknown control statement \'"<< prefix << ":" << word << "\' on line " << line << " of input file" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      return false;
   }
} // end of namespace gpu
