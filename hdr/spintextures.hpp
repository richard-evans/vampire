//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Ricardo Rama-Eiroa 2025. All rights reserved.
//
//   Email: ricardo.rama-eiroa@ed.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SPINTEXTURES_H_
#define SPINTEXTURES_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "spintextures.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for spintextures module
//--------------------------------------------------------------------------------
namespace spintextures{

   //-----------------------------------------------------------------------------
   // Function to initialise spintextures module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spintextures module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of spintextures namespace

#endif //SPINTEXTURES_H_
