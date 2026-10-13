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
#ifndef UNITS_H_
#define UNITS_H_

#include <string>
#include <vector>

namespace units{

	extern const double pi;

	// Holds all the data for one unit.
	// Keeping these three values in a struct means they always travel together —
	// there is no risk of the name, type, and conversion factor getting out of
	// step with each other (which was easy to do with the old parallel arrays).
	struct unit_t {
		std::string name;        // string the user writes in the input file, e.g. "nm"
		std::string type;        // physical quantity this unit measures, e.g. "length"
		double      conversion;  // multiply the user's value by this to get internal units
	};

	// conversion functions — public API is unchanged; all existing callers still compile
	extern int  convert(std::string input_unit, double& value, std::string& type);
	extern void convert(std::string input_unit, std::vector<double>& value, std::string& type);
	extern int  revert (std::string output_unit, double& value, std::string& type);

}


#endif /* UNITS_H_ */
