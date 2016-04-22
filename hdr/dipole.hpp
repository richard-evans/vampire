//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef DIPOLE_H_
#define DIPOLE_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "dipole.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for dipole module
//--------------------------------------------------------------------------------
namespace dipole{

   //-----------------------------------------------------------------------------
   // Function to initialise dipole module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for dipole module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of dipole namespace

#endif //DIPOLE_H_
