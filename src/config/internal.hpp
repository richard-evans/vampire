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

// config module headers
#include "internal.hpp"

namespace vout{
   // function headers
   void atoms();
   void atoms_coords();
   void cells();
   void cells_coords();
}


namespace config{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of config namespace

#endif //CONFIG_INTERNAL_H_
