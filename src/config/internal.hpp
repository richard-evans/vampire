//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) rory.pond 2016. All rights reserved.
//
//   Email: rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef CONFIG_INTERNAL_H_
#define CONFIG_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the config module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "config.hpp"

namespace config{


   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
      extern bool output_atoms_config;
      extern int output_atoms_config_rate;

      extern bool output_cells_config;
      extern int output_cells_config_rate;

      extern double field_output_min_1;
      extern double field_output_max_1;
      extern double field_output_min_2;
      extern double field_output_max_2;

      extern double atoms_output_min[3];
      extern double atoms_output_max[3];
      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void atoms();
      void atoms_coords();
      void cells();
      void cells_coords();

   } // end of internal namespace

} // end of config namespace

#endif //CONFIG_INTERNAL_H_
