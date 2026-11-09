//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) rory.pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef CONFIG_H_
#define CONFIG_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers

//--------------------------------------------------------------------------------
// Namespace for variables and functions for config module
//--------------------------------------------------------------------------------
namespace config{

   //-----------------------------------------------------------------------------
   // Function to initialise config module
   //-----------------------------------------------------------------------------
   void output();

   // true if cell configurations (cells-*.cfg, which include the cell dipole fields) are written
   bool outputs_cells();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for config module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of config namespace

#endif //CONFIG_H_
