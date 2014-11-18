//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "ltmp.hpp"

// Localised temperature pulse headers
#include "internal.hpp"

namespace ltmp{

      //-----------------------------------------------------------------------------
      // Function to check localised temperature pulse is enabled and initialised
      //-----------------------------------------------------------------------------
      bool is_enabled(){
         if(ltmp::internal::enabled && ltmp::internal::initialised) return true;
         else return false;
      }

} // end of namespace ltmp
