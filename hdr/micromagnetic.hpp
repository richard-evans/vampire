//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MICROMAGNETIC_H_
#define MICROMAGNETIC_H_

// C++ standard library headers
#include <string>
#include <vector>
// Vampire headers
#include "micromagnetic.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for micromagnetic module
//--------------------------------------------------------------------------------
namespace micromagnetic{

   extern bool discretisation_micromagnetic;
   //-----------------------------------------------------------------------------
   // Function to initialise micromagnetic module
   //-----------------------------------------------------------------------------
//   void initialize();

   extern bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   extern bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);


} // end of micromagnetic namespace

#endif //MICROMAGNETIC_H_
