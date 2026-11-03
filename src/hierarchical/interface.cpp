//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cstdlib>
#include <iostream>
#include <string>

// Vampire headers
#include "hierarchical.hpp"
#include "errors.hpp"
#include "vio.hpp"

// hierarchical module headers
#include "internal.hpp"

namespace hierarchical{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for hierarchical module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="hierarchical";
      if(key!=prefix) return false;

      //--------------------------------------------------------------------
      std::string test="output-p2p-diagnostic";
      if(word==test){
         hierarchical::internal::output_p2p_diagnostic = true;
         return true;
      }
      //--------------------------------------------------------------------
      test="moment-exchange";
      if(word==test){
         test="all-reduce";
         if(value==test){
            hierarchical::internal::moment_exchange = hierarchical::internal::moment_exchange_t::all_reduce;
            return true;
         }
         test="point-to-point";
         if(value==test){
            hierarchical::internal::moment_exchange = hierarchical::internal::moment_exchange_t::point_to_point;
            return true;
         }
         terminaltextcolor(RED);
         std::cerr << "Error: Value for \'" << prefix << ":" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"all-reduce\"" << std::endl;
         std::cerr << "\t\"point-to-point\"" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      //--------------------------------------------------------------------
      test="moment-exchange-level";
      if(word==test){
         int level = atoi(value.c_str());
         vin::check_for_valid_int(level, word, line, prefix, 0, 100, "input", "0 - 100");
         hierarchical::internal::moment_exchange_level = level;
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

} // end of hierarchical namespace
