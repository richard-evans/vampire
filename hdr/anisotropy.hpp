//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland 2016. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ANISOTROPY_H_
#define ANISOTROPY_H_

// C++ standard library headers
//#include <vector>
//#include <string>

// Vampire headers
//#include "anisotropy.hpp"
#include "material.hpp"

//--------------------------------------------------------------------------------
// namespace for variables and functions for anisotropy module
//--------------------------------------------------------------------------------
namespace anisotropy
{
   //-----------------------------------------------------------------------------
   // function to initialize anisotropy module
   //-----------------------------------------------------------------------------
   void initialize(const unsigned int num_atoms,
                   std::vector<int>& atom_type_array,
                   std::vector<double>& mu_s_array);

   void fields(std::vector<double>& spin_array_x,
               std::vector<double>& spin_array_y,
               std::vector<double>& spin_array_z,
               std::vector<int>&    type_array,
               std::vector<double>& field_array_x,
               std::vector<double>& field_array_y,
               std::vector<double>& field_array_z,
               const int start_index,
               const int end_index,
               const double temperature);

   double single_spin_energy(const int atom, const int imaterial, const double sx, const double sy, const double sz, const double temperature);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for anisotropy module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, const int max_materials);

} // end of anisotropy namespace

#endif //ANISOTROPY_H_
