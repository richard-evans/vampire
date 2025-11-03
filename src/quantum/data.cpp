//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "quantum.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside quantum module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

   } // end of internal namespace

} // end of quantum namespace

