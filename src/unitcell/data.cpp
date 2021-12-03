//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack Collings 2021. All rights reserved.
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

      int max_unit_cell_material = 100;

      // Parameters to be taken in at input for exchange function functionality      
      exchange_function_t exchange_function = nearest_neighbour;
      double exchange_interaction_range = 1.0;
      double exchange_decay = 0.4; // Angstroms
      double exchange_multiplier = 1.0;
      double exchange_shift = 0.0;
      double RKKYkf = 1.0;
      std::vector <std::vector <exchange_parameters_t> > material_exchange_parameters(max_unit_cell_material, std::vector <exchange_parameters_t>(max_unit_cell_material));
      std::vector <std::vector <double> > nn_cutoff_range(max_unit_cell_material, std::vector <double>(max_unit_cell_material, 1.0));
      std::vector <std::vector <double> > interaction_cutoff_range(max_unit_cell_material, std::vector <double>(max_unit_cell_material, 1.0));

      bool sublattice_materials = false; // flag to enable identification of atoms in simple crystals by material

   } // end of internal namespace

} // end of unitcell namespace
