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
#ifndef CELLS_INTERNAL_H_
#define CELLS_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the cells module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

// Vampire headers
#include "cells.hpp"
#include "material.hpp"

namespace cells{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
//      extern bool enabled; // enable localised temperature pulse calculation
      extern bool initialised; /// flag set if initialised
      //extern int num_cells; /// number of macro-cells

      extern std::vector<double> spin_array_x;
      extern std::vector<double> spin_array_y;
      extern std::vector<double> spin_array_z;
      extern std::vector<int> atom_type_array;
     // extern std::vector<int> atom_cell_id_array;
      extern int num_atoms;

      // shared cell arrays (shared.cpp): from the local cells of each process only (compact), or copied from the full
      // size replicated arrays
      void initialise_shared_arrays_compact(const int num_local_atoms, const std::vector<double>& atom_coords_x,
                                            const std::vector<double>& atom_coords_y, const std::vector<double>& atom_coords_z,
                                            const int total_atoms_non_filler, const double unit_cell_size_x,
                                            const double unit_cell_size_y, const double unit_cell_size_z);
      void copy_full_arrays_to_shared();
      //extern int num_local_atoms;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of cells namespace

#endif //CELLS_INTERNAL_H_
