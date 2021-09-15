//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef HIERARCHICAL_H_
#define HIERARCHICAL_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "hierarchical.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for hierarchical module
//--------------------------------------------------------------------------------
namespace hierarchical{

   //-----------------------------------------------------------------------------
   // Function to initialise hierarchical module
   //-----------------------------------------------------------------------------
   void initialize(const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   std::vector<double>& atom_coords_x, //atomic coordinates
                   std::vector<double>& atom_coords_y,
                   std::vector<double>& atom_coords_z,
                   int num_atoms);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for hierarchical module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //---------------------------------------------------------------------------
   // Function to update hierarchical dipole fields
   //---------------------------------------------------------------------------
   void update(std::vector <double>& x_spin_array, // atomic spin directions
               std::vector <double>& y_spin_array,
               std::vector <double>& z_spin_array,
               std::vector <double>& m_spin_array, // atomic spin moment
               std::vector < bool >& magnetic);    // is magnetic);

} // end of hierarchical namespace

#endif //HIERARCHICAL_H_
