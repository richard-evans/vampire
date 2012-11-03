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
/// @brief Contains the vmath header namespace for sundry math utilities in zspin. 
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    04/03/2010
/// @internal
///	Created:		04/03/2010
///	Revision:	  ---
///=====================================================================================
///

#include<vector>
#include <cmath>

/// @namespace ns
/// @brief vmath namespace containing sundry math functions for zspin.
/// 
/// @internal
///=====================================================================================
///
namespace vmath{

	/// @brief Function to determine if point is within a polygon.
	///
	/// @section License
	/// Use of this code, either in source or compiled form, is subject to license from the authors.
	/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
	///
	/// @section Information
	/// @author  Richard Evans, rfle500@york.ac.uk
	/// @version 1.0
	/// @date    04/03/2010
	///
	/// @param[in] x x coordinate of point to be tested
	/// @param[in] y y coordinate of point to be tested
	/// @param[in] *polyX array of x coordinate points for the polygon
	/// @param[in] *polyY array of y coordinate points for the polygon
	/// @param[in] nsides number of sides of the polygon
	/// @return variable returned from the function
	///
	/// @internal
	///	Created:		04/03/2010
	///	Revision:	  ---
	///=====================================================================================
	///
	extern bool point_in_polygon(double, double, double*, double*, int);
	extern bool point_in_polygon2(double, double, std::vector<double>&, std::vector<double>&, int);
	extern std::vector<std::vector<double> > matmul(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);
	extern std::vector<double> matmul(std::vector<double> &, std::vector<std::vector<double> > &);
	extern std::vector<std::vector<double> > transpose(std::vector<std::vector<double> > &);
	extern std::vector<std::vector<double> > set_matrix(const unsigned int, const unsigned int, std::vector<double> &);
	extern std::vector<std::vector<double> > set_matrix(const unsigned int, const unsigned int);
	extern void print_matrix( std::vector<std::vector<double> >&);
	extern void set_rotational_matrix(double, double, double, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &, std::vector< std::vector<double> > &);
	
	extern double sign(double);
	extern int sign(int);
	inline int iround( double value ){
		return static_cast<int>(floor( value + 0.5 ));
	}
	inline int iceil( double value ){
		return static_cast<int>(ceil( value ));
	}
	
	
	
}


