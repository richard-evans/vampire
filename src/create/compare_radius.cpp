//-----------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "create.hpp"

// Internal create header
#include "internal.hpp"

namespace create{

   namespace internal{

      // comparison function for reverse order sorting
      bool compare_radius(core_radius_t first,core_radius_t second){
         if(first.radius<second.radius) return false;
         else return true;
      }

   } // end of internal namespace

} // end of create namespace
