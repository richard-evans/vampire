//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "dipole.hpp"
#include "errors.hpp"
#include "vio.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for dipole module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="dipole";
      if(key!=prefix) return false;

      //-------------------------------------------------------------------
      std::string test="solver";
      if(word==test){
         // enable dipole calculation
         dipole::activated=true;
         return true;
      }

      test="solver-fft";
      if(word==test){
         // enable dipole calculation
         dipole::fft=true;
         return true;
      }

      test="field-update-rate";
      if(word==test){
         int dpur=atoi(value.c_str());
         vin::check_for_valid_int(dpur, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         dipole::update_rate=dpur;
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

} // end of dipole namespace
