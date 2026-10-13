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

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

//--------------------------------------------------------------------------------
// Namespace for variables and functions for exchange module
//--------------------------------------------------------------------------------
namespace constants{

   // Fundamental constants (CODATA 2018 / 2019 SI redefinition)
   extern const double muB;     // Bohr magneton              (J/T)
   extern const double kB;      // Boltzmann constant         (J/K)
   extern const double kB_eV;   // Boltzmann constant         (eV/K)
   extern const double eV;      // Elementary charge / 1 eV   (J)
   extern const double gamma_SI;// Electron gyromagnetic ratio (rad/s/T)

   // Derived constants (computed from the fundamentals above)
   extern const double i_muB;           // 1/muB — frequently needed in MC and dipole code (J/T)^-1
   extern const double dipole_prefactor;// mu_0*muB/(4*pi*Ang^3) = 1e-7 * muB / 1e-30  (T)

} // end of constants namespace

#endif //CONSTANTS_H_
