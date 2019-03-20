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

#ifndef SPINTRANSPORT_H_
#define SPINTRANSPORT_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "spintransport.hpp"

// Alias spin_transport to st namespace for code brevity
namespace spin_transport = st;

//--------------------------------------------------------------------------------
// Namespace for variables and functions for spintransport module
//--------------------------------------------------------------------------------
namespace spin_transport{

   //-----------------------------------------------------------------------------
   // Function to initialise spintransport module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spintransport module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of spin_transport namespace

#endif //SPINTRANSPORT_H_
