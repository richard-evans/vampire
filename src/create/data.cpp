//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------

      namespace internal{

         //----------------------------------------------------------------------------
         // Shared variables used within create module
         //---------------------------------------------------------------------------
         std::vector<create::internal::mp_t> mp; // array of material properties

      } // end of internal namespace

} // end of create namespace
