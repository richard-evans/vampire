//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2016. All rights reserved.
//
//   Email: am1808@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef CELLS_H_
#define CELLS_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers

//--------------------------------------------------------------------------------
// Namespace for variables and functions for cells module
//--------------------------------------------------------------------------------
namespace cells{

   extern int num_atoms_in_unit_cell;
   extern int num_cells; /// number of macro-cells
   extern int num_local_cells; /// number of macro-cells
   extern double macro_cell_size; /// lateral size of local macro-cells (A)
   extern double macro_cell_size_x; /// lateral size of local macro-cells (A)
   extern double macro_cell_size_y; /// lateral size of local macro-cells (A)
   extern double macro_cell_size_z; /// lateral size of local macro-cells (A)

   extern std::vector <int> local_cell_array;
   extern std::vector <int> num_atoms_in_cell; /// number of atoms in each cell
   extern std::vector <int> num_atoms_in_cell_global; /// global number of atoms in each cell
   extern std::vector < std::vector <int> > index_atoms_array;
   extern std::vector<int> index_atoms_array1D;

   extern std::vector<double> volume_array;
   extern std::vector<double> cell_coords_array_x; /// arrays to store cells positions
   extern std::vector<double> cell_coords_array_y;
   extern std::vector<double> cell_coords_array_z;
   extern std::vector < std::vector <double> > atom_in_cell_coords_array_x;
   extern std::vector < std::vector <double> > atom_in_cell_coords_array_y;
   extern std::vector < std::vector <double> > atom_in_cell_coords_array_z;

   extern std::vector<int> cell_id_array;
   extern std::vector<int> atom_cell_id_array;
   extern std::vector<double> mag_array_x; /// arrays to store cells magnetisation
   extern std::vector<double> mag_array_y;
   extern std::vector<double> mag_array_z;
   extern std::vector<double> field_array_x; /// arrays to store cells field
   extern std::vector<double> field_array_y;
   extern std::vector<double> field_array_z;

   extern std::vector<double> pos_and_mom_array;
   extern std::vector<double> pos_array;

   extern std::vector < double > num_macro_cells_fft; /// lateral size of local macro-cells (A)
   extern std::vector<double> fft_cell_id_array;

   //---------------------------------------------------------------------------
   // Function to calculate magnetisation in cells
   //---------------------------------------------------------------------------
   extern int mag();

   //-----------------------------------------------------------------------------
   // Function to initialise cells module
   //-----------------------------------------------------------------------------
   void initialize(const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const double unit_cell_size_x,
                   const double unit_cell_size_y,
                   const double unit_cell_size_z,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const std::vector<int>& atom_cell_id_array,
                   const int num_total_atoms_for_dipole,
                   const int num_atoms
   );

   //---------------------------------------------------------------------------
   // Function to process input file parameters for cells module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of cells namespace

#endif //CELLS_H_
