//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains functions to calculate energy for a spin/system
///
/// @details None
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation.
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///

// Standard Libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire Header files
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "dipole.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "spintransport.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// sim module header
#include "internal.hpp"

namespace sim{



/// @brief Calculates the applied field energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom
/// @param[in] Sy y-spin of local atom
/// @param[in] Sz z-spin of local atom
/// @return applied field energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
double spin_applied_field_energy(const double Sx, const double Sy, const double Sz){

	return -sim::H_applied*(sim::H_vec[0]*Sx + sim::H_vec[1]*Sy + sim::H_vec[2]*Sz);

}

double spin_local_applied_field_energy(const int mat, const double sx, const double sy, const double sz){

	const double B = mp::material[mat].applied_field_strength;
	const double hx = B * mp::material[mat].applied_field_unit_vector[0];
	const double hy = B * mp::material[mat].applied_field_unit_vector[1];
	const double hz = B * mp::material[mat].applied_field_unit_vector[2];

	return -(hx*sx + hy*sy + hz*sz);

}

/// @brief Calculates the magnetostatic energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2012. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    08/06/2012
///
/// @param[in] atom atom number
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom
/// @param[in] Sy y-spin of local atom
/// @param[in] Sz z-spin of local atom
/// @return magnetostatic energy
///
/// @internal
///	Created:		08/06/2012
///	Revision:	  ---
///=====================================================================================
///
double spin_magnetostatic_energy(const int atom, const double Sx, const double Sy, const double Sz){
   //return -1.0*(dipole::atom_dipolar_field_array_x[atom]*Sx+dipole::atom_dipolar_field_array_y[atom]*Sy+dipole::atom_dipolar_field_array_z[atom]*Sz);
   return -1.0*(dipole::atom_mu0demag_field_array_x[atom]*Sx+dipole::atom_mu0demag_field_array_y[atom]*Sy+dipole::atom_mu0demag_field_array_z[atom]*Sz);
}

/// @brief Calculates the total energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number
/// @return total spin energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
double calculate_spin_energy(const int atom){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "calculate_spin_energy has been called" << std::endl;

	// Local spin value
	const double Sx=atoms::x_spin_array[atom];
	const double Sy=atoms::y_spin_array[atom];
	const double Sz=atoms::z_spin_array[atom];

	// Determine neighbour material
	const int imaterial=atoms::type_array[atom];

	// Initialise energy to zero
	double energy=0.0;

	// Calculate total spin energy
   energy += exchange::single_spin_energy(atom, Sx, Sy, Sz);
   energy += exchange::single_spin_biquadratic_energy(atom, Sx, Sy, Sz);

   // calculate anisotropy energy for atom
   energy += anisotropy::single_spin_energy(atom, imaterial, Sx, Sy, Sz, sim::temperature);

	energy+=spin_applied_field_energy(Sx, Sy, Sz);
	energy+=spin_magnetostatic_energy(atom, Sx, Sy, Sz);

	// local applied fields
	if(sim::local_applied_field) energy += spin_local_applied_field_energy(imaterial, Sx, Sy, Sz);

	// vcma energy
	const double vcma = program::fractional_electric_field_strength * spin_transport::get_voltage() * sim::internal::vcmak[imaterial];
	energy -= vcma * Sz * Sz;

	return energy; // Tesla
}

} // end of namespace sim
