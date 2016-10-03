//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

//   bool micro;

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside micromagnetic module
      //------------------------------------------------------------------------

      std::vector<double> ms;
      std::vector<double> ku_z;
      std::vector<double> A;
      std::vector<double> Chi;
      std::vector<double> tc;
      //bool discretisation_micromagnetic(false);

   } // end of internal namespace

} // end of micromagnetic namespace
