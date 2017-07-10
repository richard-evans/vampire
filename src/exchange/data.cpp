//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside exchange module
      //------------------------------------------------------------------------
      std::vector<internal::mp_t> mp; // array of material properties

      bool enable_dmi = false; // flag to enable dmi calculation

      double dmi_cutoff_range = 2.6; // cutoff range for DMI calculation (Ångstroms)

      exchange_t exchange_type = isotropic; // exchange type to use in simulation

      bool use_material_exchange_constants = true; // flag to enable material exchange parameters

   } // end of internal namespace

} // end of exchange namespace
