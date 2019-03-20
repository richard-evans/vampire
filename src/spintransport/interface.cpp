//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "spintransport.hpp"
#include "errors.hpp"
#include "vio.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spintransport module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-transport";
      if(key!=prefix) return false;

      // cell-size
      // current-direction
      // channel length
      // voltage

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

      // tunnel-barrier = true/false (if sp, then propagate across)
      // spin-resistivity (Ohm/m^2) [AP state] - tunnel barrier has a high R and high SR, normal materials SR is low ~ a few %
      // resistivity (Ohm/m^2) [P state]
      // specific-heat-capacity? for joule heating

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of spin_transport namespace
