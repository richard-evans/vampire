//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spintransport.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside spintransport module
      //------------------------------------------------------------------------

      // enumerated list of different current directions
      enum current_direction_t {px,py,pz,mx,my,mz}; // +x,+y,+z,-x,-y,-z
      current_direction_t current_direction = pz; // current direction (default along +z direction)

      double cell_size_x = 10.0; // cell size along x-direction (1 nm default size)
      double cell_size_y = 10.0; // cell size along y-direction (1 nm default size)
      double cell_size_z = 10.0; // cell size along z-direction (1 nm default size)

      // array of stacks (list of cells) along current direction
      std::vector < std::vector <unsigned int> > stack_array;

   } // end of internal namespace

} // end of spin_transport namespace
