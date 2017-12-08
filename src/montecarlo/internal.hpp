//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MONTECARLO_INTERNAL_H_
#define MONTECARLO_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the montecarlo module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "montecarlo.hpp"

// montecarlo module headers
#include "internal.hpp"

namespace montecarlo{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      //MC Variables
      extern std::vector<double> Sold;
      extern std::vector<double> Snew;

      //MC-MPI variables
      extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant
      extern bool mc_parallel_initialized;
      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void mc_move(const std::vector<double>&, std::vector<double>&);

   } // end of internal namespace

} // end of montecarlo namespace

#endif //MONTECARLO_INTERNAL_H_
