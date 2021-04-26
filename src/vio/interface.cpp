//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include<algorithm> 

// Vampire headers
#include "errors.hpp"
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

namespace vio{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for vio module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="output";
      if(key!=prefix) return false;

      //-------------------------------------------------------------------
      std::string test="precision";
      if(word==test){
         int p=atoi(value.c_str());
         vin::check_for_valid_int(p, word, line, prefix, 4, 16,"input","4 - 16");
         vout::custom_precision = true; // enable user definable precision for output
         vout::precision = p;
         return true;
      }
      //-------------------------------------------------------------------
      test="fixed-width";
      if(word==test){
         vout::fixed=true;
         vout::fw_size = std::max(int(vout::precision+7),vout::max_header); // enable fixed width output
         vout::fw_size_int = std::max(int(vout::precision),vout::max_header); // enable fixed width output
         return true;
      }
      test="column-headers";
      if(word==test){
        vout::header_option = true;
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

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of vio namespace
