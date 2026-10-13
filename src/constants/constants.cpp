//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers (need this for definition of external linkage of const variables)
#include "constants.hpp"

//--------------------------------------------------------------------------------
// Namespace for program constants
//
// All values are taken from CODATA 2018 / the 2019 SI redefinition.
// eV, kB, and muB are now exact by definition under the 2019 SI.
// This is the single source of truth for physical constants in vampire —
// do not hard-code these values elsewhere in the codebase.
//--------------------------------------------------------------------------------
namespace constants{

   // Fundamental constants (CODATA 2018)
   const double muB      = 9.2740100783e-24; // Bohr magneton              (J/T)
   const double kB       = 1.380649e-23;     // Boltzmann constant         (J/K)   exact
   const double kB_eV    = 8.617333262e-5;   // Boltzmann constant         (eV/K)
   const double eV       = 1.602176634e-19;  // 1 electron volt            (J)     exact
   const double gamma_SI = 1.76085963023e11; // Electron gyromagnetic ratio (rad/s/T)

   // Derived constants — computed here so the relationship to the fundamentals
   // is explicit and the values stay consistent if the fundamentals are updated.
   const double i_muB            = 1.0 / muB;               // (J/T)^-1
   const double dipole_prefactor = 1.0e-7 * muB / 1.0e-30; // mu_0*muB/(4*pi*Ang^3) (T)

} // end of constants namespace
