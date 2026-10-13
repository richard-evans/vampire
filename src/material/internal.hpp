//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MATERIAL_INTERNAL_H_
#define MATERIAL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the material module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>

// Vampire headers
#include "material.hpp"

namespace mp{

   namespace internal{

      //-------------------------------------------------------------------------
      // Function protypes for functions inside: read.cpp
      //-------------------------------------------------------------------------
      int read_mat_file(std::string const, int const);

      //-------------------------------------------------------------------------
      // Function protypes for functions inside: match.cpp
      //-------------------------------------------------------------------------
      int match_material(std::string const, std::string const, std::string const, int const, int const, int const, std::string const, std::string const);

   } // end of internal namespace

} // end of mp namespace

#endif //MATERIAL_INTERNAL_H_
