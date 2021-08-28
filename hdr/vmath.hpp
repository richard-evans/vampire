//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#include <vector>
#include <cmath>

/// @namespace ns
/// @brief vmath namespace containing sundry math functions for vampire.
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
	extern bool point_in_polygon_factor(double, double, double, double*, double*, int);
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
   // rounding function for 64 bit integers
   inline int64_t iround64( double value ){
		return static_cast<int64_t>(floor( value + 0.5 ));
	}
	inline int iceil( double value ){
		return static_cast<int>(ceil( value ));
	}

	extern double interpolate_m(double,double,double,double);
	extern double interpolate_c(double,double,double,double);
	extern double minimum3(double first, double second, double third);
	extern void regression(std::vector<double>& x, std::vector<double>& y, double& m, double& c);

}
