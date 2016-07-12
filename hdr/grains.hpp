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
#ifndef GRAINS_H_
#define GRAINS_H_

#include <vector>
#include <fstream>

/// @namespace
/// @brief Contains data about all grains in the system.
///
/// @internal
///=====================================================================================
///
namespace grains{

	extern int num_grains;
	extern bool random_anisotropy; // flag to control randomly oriented uniaxial anisotropy

	extern std::vector <int> grain_size_array;

	extern std::vector <double> x_coord_array;
	extern std::vector <double> y_coord_array;
	extern std::vector <double> z_coord_array;

	extern std::vector <double> x_mag_array;
	extern std::vector <double> y_mag_array;
	extern std::vector <double> z_mag_array;
	extern std::vector <double> mag_m_array;

	extern std::vector <double> sat_mag_array;

	extern int set_properties();
	extern int mag();
	extern int output_mag(std::ofstream&);
	extern void output_mat_mag(std::ostream&);
}

#endif /*GRAINS_H_*/
