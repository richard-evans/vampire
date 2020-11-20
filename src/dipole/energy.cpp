//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2018. All rights reserved.
//
//------------------------------------------------------------------------------
//

// dipole module headers
#include "internal.hpp"

namespace dipole{

//------------------------------------------------------------------------------
// simple function to calculate energy of spin in dipole (magnetostatic) field
//------------------------------------------------------------------------------
double spin_magnetostatic_energy(const int atom, const double sx, const double sy, const double sz){

   return -1.0 * ( dipole::atom_mu0demag_field_array_x[atom] * sx + dipole::atom_mu0demag_field_array_y[atom] * sy + dipole::atom_mu0demag_field_array_z[atom] * sz);

}

} // end of dipole namespace
