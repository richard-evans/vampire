//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "quantum.hpp"
#include "errors.hpp"
#include "vio.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for quantum module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="quantum";
      if(key!=prefix) return false;

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      std::string test = "noise-type";
      if( word == test ){
            test = "classical";
            if(value == test){
               noise_type = classical;
               return true;
            }
            test = "quantum";
            if( value == test ){
               noise_type = quantum_zero;
               return true;
            }
            test = "quantum-no-zero";
            if( value == test ){
               noise_type = quantum_no_zero;
               return true;
            }
            else{
               terminaltextcolor(RED);
               std::cerr << "Error - value for \'quantum:" << word << "\' must be one of:" << std::endl;
               std::cerr << "\t\"classical\"" << std::endl;
               std::cerr << "\t\"quantum\"" << std::endl;
               std::cerr << "\t\"quantum-no-zero\"" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }
         }
      }
      //------------------------------------------------------------------------
      test = "noise-window-size";
      if( word == test ){
         int ws = vin::str_to_int(value);
         vin::check_for_valid_int(ws, word, line, prefix, 1, 6000000000, "input", "> 0");
         if(ws % 6 != 0){
            std::cerr << "Error: Quantum window size must be divisible by 6, eg 6000, 60000 etc." << std::endl;
            return false;
         }
         internal::window_size = ws;
         return true;
      }
      //------------------------------------------------------------------------
      test = "noise-decimation-factor";
      if( word == test ){
         int md = vin::str_to_int(value);
         vin::check_for_valid_int(md, word, line, prefix, 1, 1000000, "input", "> 0");
         internal::M_decimation = md;
         return true;
      }
      //------------------------------------------------------------------------

      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Resize mp vector if necessary
      if(internal::mp.size() <= (unsigned int)super_index){
         internal::mp.resize(super_index + 1);
      }

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      std::string test = "quantum-A";
      if( word == test ){
         double A = vin::str_to_double(value);
         vin::check_for_valid_value(A, word, line, prefix, unit, "energy", 0.0, 1.0e9, "material", "> 0");
         internal::mp[super_index].A.set(A);
         return true;
      }

      test = "quantum-Gamma";
      if( word == test ){
         double Gamma = vin::str_to_double(value);
         vin::check_for_valid_value(Gamma, word, line, prefix, unit, "frequency", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].Gamma.set(Gamma);
         return true;
      }

      test = "quantum-omega0";
      if( word == test ){
         double omega0 = vin::str_to_double(value);
         vin::check_for_valid_value(omega0, word, line, prefix, unit, "frequency", 0.0, 1.0e15, "material", "> 0");
         internal::mp[super_index].omega0.set(omega0);
         return true;
      }

      test = "quantum-S0";
      if( word == test ){
         double S0 = vin::str_to_double(value);
         vin::check_for_valid_value(S0, word, line, prefix, unit, "none", 0.0, 100.0, "material", "> 0");
         internal::mp[super_index].S0.set(S0);
         return true;
      }

      return false;

   }

} // end of quantum namespace
