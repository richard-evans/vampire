//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "montecarlo.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for montecarlo module
//--------------------------------------------------------------------------------
namespace montecarlo{

   //-----------------------------------------------------------------------------
   // Function to initialise montecarlo module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for montecarlo module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of montecarlo namespace

#endif //MONTECARLO_H_
