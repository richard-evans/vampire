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
#ifndef CELLS_H_
#define CELLS_H_

#include <vector>
#include <fstream>

/// @namespace
/// @brief Contains data about all cells in the system.
///
/// @internal
///=====================================================================================
///
namespace cells{

	extern int num_cells;
	extern int num_local_cells;
   extern int num_atoms_in_unit_cell;

	extern double size;

	extern bool initialised;

	extern std::vector <int> num_atoms_in_cell;
	extern std::vector <int> local_cell_array;

	extern std::vector <double> x_coord_array;
	extern std::vector <double> y_coord_array;
	extern std::vector <double> z_coord_array;

	extern std::vector <double> x_mag_array;
	extern std::vector <double> y_mag_array;
	extern std::vector <double> z_mag_array;
	
	extern std::vector <double> x_field_array;
	extern std::vector <double> y_field_array;
	extern std::vector <double> z_field_array;

   extern std::vector <double> volume_array;

	extern int initialise();
	extern int mag();
	extern int output_mag(std::ofstream&);
}

#endif /*CELLS_H_*/


