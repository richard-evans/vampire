//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Ricardo Rama-Eiroa 2025. All rights reserved.
//
//   Email: ricardo.rama-eiroa@ed.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spintextures.hpp"

// spintextures module headers
#include "internal.hpp"

namespace spintextures{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside spintextures module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

   } // end of internal namespace

} // end of spintextures namespace

