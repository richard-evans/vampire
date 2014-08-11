//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "errors.hpp"
#include "ltmp.hpp"
#include "vio.hpp"

// Localised temperature pulse headers
#include "internal.hpp"

namespace ltmp{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ltmp settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="local-temperature-pulse";
      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="cell-size";
      if(word==test){
         double csize=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(csize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         ltmp::internal::micro_cell_size = csize;
         return true;
      }
      //--------------------------------------------------------------------
      test="laser-spot-size";
      if(word==test){
         double lssize=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(lssize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         ltmp::internal::laser_spot_size = lssize;
         return true;
      }
      //--------------------------------------------------------------------
      test="penetration-depth";
      if(word==test){
         double pdepth=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(pdepth, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         ltmp::internal::penetration_depth = pdepth;
         return true;
      }
      //--------------------------------------------------------------------
      test="thermal-conductivity";
      if(word==test){
         double tdc=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(tdc, word, line, prefix, unit, "none", 0.0, 100,"input","0.0 - 100.0 J/m/s/K");
         ltmp::internal::thermal_conductivity = tdc;
         return true;
      }
      //--------------------------------------------------------------------
      test="output-microcell-data";
      if(word==test){
         ltmp::internal::output_microcell_data = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="temperature-profile";
      if(word==test){
         test="lateral";
         if(value==test){
            ltmp::internal::lateral_discretisation = true;
            ltmp::internal::vertical_discretisation = false;
            ltmp::internal::enabled = true;
            return true;
         }
         test="vertical";
         if(value==test){
            ltmp::internal::lateral_discretisation = false;
            ltmp::internal::vertical_discretisation = true;
            ltmp::internal::enabled = true;
            return true;
         }
         test="lateral-vertical";
         if(value==test){
            ltmp::internal::lateral_discretisation = true;
            ltmp::internal::vertical_discretisation = true;
            ltmp::internal::enabled = true;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"lateral\"" << std::endl;
            std::cerr << "\t\"vertical\"" << std::endl;
            std::cerr << "\t\"lateral-vertical\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
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
} // end of namespace ltmp

//material[1]:kerr-sensitivity-depth (fit to exp (-z/factor) )
