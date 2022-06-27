//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SLD_H_
#define SLD_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "sld.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for sld module
//--------------------------------------------------------------------------------
namespace sld{

   //-----------------------------------------------------------------------------
   // Function to initialise sld module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for sld module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of sld namespace

#endif //SLD_H_
