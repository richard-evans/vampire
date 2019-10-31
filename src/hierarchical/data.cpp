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

// C++ standard library headers

// Vampire headers
#include "hierarchical.hpp"

// hierarchical module headers
#include "internal.hpp"
#include <vector>

namespace hierarchical{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside hierarchical module
      //------------------------------------------------------------------------

      int num_levels; // number of levels for hierarchical expansion

      //create arrays for data storage
      std::vector < double > cell_positions;             // positions of hierarchical cells 3N
      std::vector < double > cell_positions_mom;         // positions and magnetic moments of hierarchical cells 4N
      std::vector < double > cell_dimensions;            // shape and size of each cell
      std::vector < int > cells_level_start_index;       // first cell in Level
      std::vector < int > cells_level_end_index;         // last cell in Level
      std::vector < int > cells_in_cells;                // 1D list of which cells in level L are contained within level L+1
      std::vector < int > interaction_range;             // a list of interaction ranges for each hierarchical level
      std::vector < int > cells_in_cells_start_index;    // index of first cell in level L-1
      std::vector < int > cells_in_cells_end_index;      // index of last  cell in level L-1
      std::vector <int> interaction_list;                // list of interacting hierarchical cells
      std::vector < int > interaction_list_start_index;  // index of first neighbour in interactions list
      std::vector < int > interaction_list_end_index;    // index of last  neighbour in interactions list
      std::vector < int > num_atoms_in_cell;             // list of number of atoms within each cell
      std::vector < int > proc_cell_index_array1D;       //
      std::vector < int > cells_num_atoms_in_cell;       //

      double av_cell_size;                               // weighted average of cell_size_x, cell_size_y, cell_size_z

      std::vector<double> mag_array_x;       // arrays to store cell magnetisation
      std::vector<double> mag_array_y;
      std::vector<double> mag_array_z;

      int total_num_cells;
      int num_zero_level_cells;              // number of zero level cells

      std::vector<double> rij_tensor_xx;     // dipole-dpole tensor coponnets for hierarchical method
      std::vector<double> rij_tensor_xy;
      std::vector<double> rij_tensor_xz;

      std::vector<double> rij_tensor_yy;
      std::vector<double> rij_tensor_yz;
      std::vector<double> rij_tensor_zz;

      std::vector<double> ha_cells_pos_and_mom_array;    //
      std::vector < int > ha_proc_cell_index_array1D;    //

   } // end of internal namespace

} // end of hierarchical namespace
