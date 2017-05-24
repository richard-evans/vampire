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
#include "sim.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------

   namespace internal{

      //----------------------------------------------------------------------------
      // Shared variables used within sim module
      //---------------------------------------------------------------------------
      std::vector<sim::internal::mp_t> mp; // array of material properties
      std::vector<double> slonczewski_aj; // array of adiabatic spin torques
      std::vector<double> slonczewski_bj; // array of non-adiabatic spin torques
      std::vector<double> slonczewski_spin_polarization_unit_vector(3,0.0); // spin polarization direction

      int num_monte_carlo_preconditioning_steps = 0;

   } // end of internal namespace

} // end of sim namespace
