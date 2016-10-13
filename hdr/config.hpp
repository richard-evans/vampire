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

namespace vout{
	extern void config();

	extern double field_output_min_1;
	extern double field_output_max_1;
	extern double field_output_min_2;
	extern double field_output_max_2;

	extern bool output_cells_config;
	extern int output_cells_config_rate;

	extern bool output_atoms_config;
	extern int output_atoms_config_rate;

	extern double atoms_output_min[3];
	extern double atoms_output_max[3];
}
//--------------------------------------------------------------------------------
// Namespace for variables and functions for config module
//--------------------------------------------------------------------------------
namespace config{

   //-----------------------------------------------------------------------------
   // Function to initialise config module
   //-----------------------------------------------------------------------------
   void initialize();

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
