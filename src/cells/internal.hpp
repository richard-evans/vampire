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

// cells module headers
#include "internal.hpp"

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

      // cells:probe-non-magnetic-cells - opt-in flag (default off) to expose
      // macrocells containing only non-magnetic atoms as valid dipole-field
      // evaluation points (e.g. a dense non-magnetic sensor/read-back layer).
      // When off, behaviour is bit-for-bit identical to before this option
      // existed.
      extern bool probe_non_magnetic_cells;

      // Number of atoms of any kind (magnetic + non-magnetic) in each cell.
      // Unlike cells::num_atoms_in_cell (which counts magnetic atoms only and
      // must keep that meaning everywhere else), this is only used to build
      // cells::cell_id_array / cells::local_cell_array / cells::pos_and_mom_array
      // for non-magnetic-only cells when probe_non_magnetic_cells is enabled.
      extern std::vector<int> num_any_atoms_in_cell;

      extern std::vector<double> total_moment_array;
      extern std::vector<double> cell_position_array; /// position of cells in x,y,z (3*n) MIRRORED on all CPUs // dont need this

      extern std::vector<double> spin_array_x;
      extern std::vector<double> spin_array_y;
      extern std::vector<double> spin_array_z;
      extern std::vector<int> atom_type_array;
     // extern std::vector<int> atom_cell_id_array;
      extern int num_atoms;
      //extern int num_local_atoms;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of cells namespace

#endif //CELLS_INTERNAL_H_
