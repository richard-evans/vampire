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

// C++ standard library headers
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>

// Vampire headers
#include "cells.hpp"

// cells module headers
#include "internal.hpp"

namespace cells{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   int num_atoms_in_unit_cell=0;
   int num_cells; /// number of macro-cells
   int num_local_cells=0; /// number of macro-cells
   double macro_cell_size = 10.0; /// macro-cells size (A)
   double macro_cell_size_x = 10.0; /// macro-cells size (A)
   double macro_cell_size_y = 10.0; /// macro-cells size (A)
   double macro_cell_size_z = 10.0; /// macro-cells size (A)


   std::vector <int> local_cell_array;
   std::vector<int> num_atoms_in_cell; /// number of atoms in each cell
   std::vector<int> num_atoms_in_cell_global; /// global number of atoms in each cell
   std::vector < std::vector <int> > index_atoms_array; /// array to store list of atoms associated with cells in 2D structure
   std::vector<int> index_atoms_array1D; /// array to store list of atoms associated with cells in 1D structure

   std::vector<double> volume_array;
   std::vector < std::vector <double> > atom_in_cell_coords_array_x;
   std::vector < std::vector <double> > atom_in_cell_coords_array_y;
   std::vector < std::vector <double> > atom_in_cell_coords_array_z;
   std::vector<int> cell_id_array;        /// array to store index of local cells
   std::vector<int> atom_cell_id_array;   /// array to store index of cells associated with each atom
   std::vector<double> mag_array_x; /// arrays to store cells magnetisation
   std::vector<double> mag_array_y;
   std::vector<double> mag_array_z;
   std::vector<double> field_array_x; /// arrays to store cells field
   std::vector<double> field_array_y;
   std::vector<double> field_array_z;

   std::vector<double> pos_and_mom_array; /// arrays to store cells positions
   std::vector<double> pos_array; /// arrays to store cells positions

   std::vector < double > num_macro_cells_fft(3,10.0); /// macro-cells size (A)
   std::vector<double> fft_cell_id_array; /// arrays to store cells positions

   //---------------------------------------------------------------------------
   // Function to calculate magnetisation in cells
   //---------------------------------------------------------------------------
   int mag();

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside cells module
      //------------------------------------------------------------------------
      //bool enabled=false; // enable localised temperature pulse calculation
      bool initialised=false; /// flag set if initialised
      std::vector<double> volume_array;
      std::vector<double> total_moment_array;
      std::vector<double> cell_position_array;

      std::vector<double> spin_array_x;
      std::vector<double> spin_array_y;
      std::vector<double> spin_array_z;
      std::vector<int> atom_type_array;
      int num_atoms;
   } // end of internal namespace

} // end of cells namespace
