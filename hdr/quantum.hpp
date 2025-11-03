//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

#ifndef QUANTUM_H_
#define QUANTUM_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "quantum.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for quantum module
//--------------------------------------------------------------------------------
namespace quantum{

   //-----------------------------------------------------------------------------
   // Function to initialise quantum module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for quantum module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of quantum namespace

#endif //QUANTUM_H_
