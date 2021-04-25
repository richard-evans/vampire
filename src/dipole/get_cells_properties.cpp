//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo Meo 2021. All rights reserved.
//
//-----------------------------------------------------------------------------


// C++ standard library headers

// Vampire headers
#include "dipole.hpp"

// exchange module headers
#include "internal.hpp"

namespace dipole{

   //------------------------------------------------------------------------------
   // Function to return internal num_atoms_in_cell vector
   //------------------------------------------------------------------------------
   std::vector<int> get_num_atoms_in_cell_array(){
      return internal::cells_num_atoms_in_cell;
   }

   //------------------------------------------------------------------------------
   // Function to return number of atoms in cell
   //------------------------------------------------------------------------------
   unsigned int get_num_atoms_in_cell(const int cell){
      return internal::cells_num_atoms_in_cell[cell];
   }

   //------------------------------------------------------------------------------
   // Function to return total number of cells
   //------------------------------------------------------------------------------
   unsigned int get_tot_num_cells(){
      return internal::cells_num_cells;
   }

   //------------------------------------------------------------------------------
   // Function to return total number of cells
   //------------------------------------------------------------------------------
   unsigned int get_tot_num_local_cells(){
      return internal::cells_num_local_cells;
   }

} // end of dipole namespace

