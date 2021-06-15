//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside unitcell module
      //------------------------------------------------------------------------
      std::string crystal_structure="sc";
      std::string unit_cell_filename="";

      double unit_cell_size_x = 3.54;
      double unit_cell_size_y = 3.54;
      double unit_cell_size_z = 3.54;

      exchange_function_t exchange_function = nearest_neighbour;
      double exchange_interaction_range = 1.0;
      double exchange_decay = 0.4; // Angstroms

      double exchange_parameter_A = 1.0;
      double exchange_parameter_B = 1.0;
      double exchange_parameter_C = 0.0;

      bool sublattice_materials = false; // flag to enable identification of atoms in simple crystals by material

   } // end of internal namespace

} // end of unitcell namespace
