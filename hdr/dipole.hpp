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
   void initialize(const int cells_num_atoms_in_unit_cell,
                   const int cells_num_cells, /// number of macrocells
                   const int cells_num_local_cells, /// number of local macrocells
                   const double cells_macro_cell_size
                   const std::vector <int>& cells_local_cell_array,
                   const std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                   const std::vector < std::vector <int> >& cells_index_atoms_array,
                   const std::vector<double>& cells_volume_array,
                   const std::vector<double>& cells_cell_coords_array_x, /// arrays to store cells positions
                   const std::vector<double>& cells_cell_coords_array_y,
                   const std::vector<double>& cells_cell_coords_array_z,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                   const std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                   const std::vector<double>& cells_mag_array_x, /// arrays to store cells magnetisation
                   const std::vector<double>& cells_mag_array_y,
                   const std::vector<double>& cells_mag_array_z,
                   const std::vector<double>& cells_field_array_x, /// arrays to store cells field
                   const std::vector<double>& cells_field_array_y,
                   const std::vector<double>& cells_field_array_z,
                   const std::vector<int>& atom_type_array,
                   const std::vector<int>& atom_cell_array,
                   const int num_atoms,
                   const std::vector<double>& atom_dipolar_field_array_x, /// arrays to store atoms dipolar field
                   const std::vector<double>& atom_dipolar_field_array_y,
                   const std::vector<double>& atom_dipolar_field_array_z
   );

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
