//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2016. All rights reserved.
//
//   Email: am1808@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "cells.hpp"
#include "errors.hpp"
#include "vio.hpp"

// cells module headers
#include "internal.hpp"

namespace cells{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for cells module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){
  // int match_dimension(std::string const word,sttd::string const value,std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="cells";
      if(key!=prefix) return false;

/*      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="macro-cell-size";
      if(word==test){
         double csize=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(csize, word, line, prefix, unit, "length", 0.1, 1.0e7,"input","0.1 Angstroms - 1 millimetre");
         cells::macro_cell_size = csize;
         return true;
      } */

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

} // end of cells namespace

