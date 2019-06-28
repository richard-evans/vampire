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

      int num_levels;
      //create arrays for data storage
      std::vector < double > cell_positions;
      std::vector < double > cell_positions_mom;
      std::vector < double > cell_dimensions;
      std::vector < int > cells_level_start_index;
      std::vector < int > cells_level_end_index;
      std::vector < int > cells_in_cells;
      std::vector < int > interaction_range;
      std::vector < int > cells_in_cells_start_index;
      std::vector < int > cells_in_cells_end_index;
      std::vector <int> interaction_list;
      std::vector < int > interaction_list_start_index;
      std::vector < int > interaction_list_end_index;
      std::vector < int > num_atoms_in_cell;
      std::vector < int > proc_cell_index_array1D;
      std::vector < int > cells_num_atoms_in_cell;



      std::vector<double> mag_array_x; /// arrays to store cells magnetisation
      std::vector<double> mag_array_y;
      std::vector<double> mag_array_z;

      int total_num_cells;
      int num_zero_level_cells;

      std::vector<double> rij_tensor_xx;
      std::vector<double> rij_tensor_xy;
      std::vector<double> rij_tensor_xz;

      std::vector<double> rij_tensor_yy;
      std::vector<double> rij_tensor_yz;
      std::vector<double> rij_tensor_zz;


      std::vector<double> ha_cells_pos_and_mom_array;
      std::vector < int > ha_proc_cell_index_array1D;

   } // end of internal namespace

} // end of hierarchical namespace
