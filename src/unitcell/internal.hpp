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

#ifndef UNITCELL_INTERNAL_H_
#define UNITCELL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the unitcell module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
	   extern std::string crystal_structure;
	   extern std::string unit_cell_filename;

      extern double unit_cell_size_x;
      extern double unit_cell_size_y;
      extern double unit_cell_size_z;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void build_simple_cubic(unitcell::unit_cell_t& unit_cell);
      void build_body_centred_cubic(unitcell::unit_cell_t& unit_cell);
      void build_face_centred_cubic(unitcell::unit_cell_t& unit_cell);
      void build_hexagonal_close_packed(unitcell::unit_cell_t& unit_cell);
      void calculate_interactions(unit_cell_t& unit_cell);
      void read_unit_cell(unit_cell_t & unit_cell, std::string filename);
      void verify_exchange_interactions(unit_cell_t & unit_cell, std::string filename);

   } // end of internal namespace

} // end of unitcell namespace

#endif //UNITCELL_INTERNAL_H_
