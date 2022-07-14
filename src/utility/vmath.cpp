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
/// @brief Contains the vmath source namespace for sundry math utilities in vampire.
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
#include "errors.hpp"
#include "vmath.hpp"
#include "vio.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace vmath{

// Following functions are adapted from http://alienryderflex.com/polygon/
//
// (c) 1998,2006,2007 Darel Rex Finley
//


bool point_in_polygon(double x, double y,double *polyX, double *polyY,int polySides) {
	///========================================================================================================
	///		 						Function to decide if point is within polygon
	///
	///														Version 1.0
	///
	///												R F Evans 15/07/2009
	///========================================================================================================

	int   i, j=polySides-1 ;
	bool  oddNodes=0;

	for (i=0; i<polySides; i++) {
		if ((polyY[i]<y && polyY[j]>=y) || (polyY[j]<y && polyY[i]>=y)) {
			if (polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x) {
				oddNodes=!oddNodes;
			}
		}
		j=i;
	}

  return oddNodes;

}

//---------------------------------------------------------------------------------------------
//
//	   						Function to determine if point is in scaled polygon
//
//												(c) R F L Evans 2015
//
//---------------------------------------------------------------------------------------------
bool point_in_polygon_factor(double x, double y, double factor, double *polyX, double *polyY,int polySides) {
	int   i, j=polySides-1 ;
	bool  oddNodes=0;

	for (i=0; i<polySides; i++) {
		if ((polyY[i]*factor<y && polyY[j]*factor>=y) || (polyY[j]*factor<y && polyY[i]*factor>=y)) {
			if (polyX[i]*factor+(y-polyY[i]*factor)/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])*factor<x) {
				oddNodes=!oddNodes;
			}
		}
		j=i;
	}

  return oddNodes;

}

bool point_in_polygon2(double x, double y, std::vector<double>& polyX, std::vector<double>& polyY, const int polySides) {
	///========================================================================================================
	///		 						Function to decide if point is within polygon
	///														Version 2.0
	///												R F Evans 11/09/2012
	///========================================================================================================

	/*// Check for correct calling of function
	bool error = false;
	if(polyX.size()!=polySides) error=true;
	if(polyY.size()!=polySides) error=true;

	// check for error
	if(error){
		zlog << zTs() << "Internal error in determining point in polygon. Number of sides " << polySides <<
		" must be equal to length of vectors " << polyX.size() << " and " << polyY.size() << ". Exiting." << std::endl;
		err::vexit();
	}*/
	x+=1e-10; // Added tiny amount to atomic coordinates to include atoms at 0,0,0
	y+=1e-10;

	int j=polySides-1 ;
	bool oddNodes=false;

	for (int i=0; i<polySides; i++) {
		if ((polyY[i]<y && polyY[j]>=y) || (polyY[j]<y && polyY[i]>=y)) {
			if (polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x) {
				oddNodes=!oddNodes;
			}
		}
		j=i;
	}

  return oddNodes;

}

/// @brief Matrix Multiplication Function
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return Matrix C
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///	Reference from Matmul guide at http://www.purplemath.com/modules/mtrxmult3.htm
///=====================================================================================
///
std::vector<std::vector<double> > matmul(
	std::vector<std::vector<double> >& A,
	std::vector<std::vector<double> >& B)
{

	// Get sizes of matrices A and B
	const unsigned int A_num_rows = A.size();
	const unsigned int A_num_cols = A[0].size();

	const unsigned int B_num_rows = B.size();
	const unsigned int B_num_cols = B[0].size();

	// Check for valid multiplication
	if(A_num_cols!=B_num_rows){
		terminaltextcolor(RED);
		std::cerr << "Error in matrix multiplication - matrices do not produce a valid product!" << std::endl;
		terminaltextcolor(WHITE);
		err::vexit();
	}

	// Set number of internal products
	const unsigned int product = A_num_cols;

	// Set sizes of result matrix C
	const unsigned int C_num_rows = A_num_rows;
	const unsigned int C_num_cols = B_num_cols;

	// Declare result matrix
	std::vector<std::vector<double> > C;
	C.resize(C_num_rows);
	for(unsigned int row=0;row<C_num_rows;row++) C[row].resize(C_num_cols,0.0);

	//std::cout << "Product is:" << std::endl;
	//C[i][j] = Ari . Bcj
	// Calculate product
	for(unsigned int i=0;i<C_num_rows;i++){
		//std::cout << "[";
		for(unsigned int j=0;j<C_num_cols;j++){
			for(unsigned int k=0; k<product ;k++){
				C[i][j]+=A[i][k]*B[k][j];
			}
			//std::cout << C[i][j] << "\t";
		}
		//std::cout << "]"<< std::endl;
	}

	// Return
	return C;
}

/// @brief Overloaded Matrix Multiplication Function for Vectors
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return Matrix C
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///	Reference from Matmul guide at http://www.purplemath.com/modules/mtrxmult3.htm
///=====================================================================================
///
std::vector<double> matmul(
	std::vector<double>& V,
	std::vector<std::vector<double> >& M)
{

	// Get sizes of vector V and matrix M
	const unsigned int V_num_cols = V.size();

	const unsigned int M_num_rows = M.size();
	const unsigned int M_num_cols = M[0].size();

	// Check for valid multiplication
	if(V_num_cols!=M_num_rows){
		terminaltextcolor(RED);
		std::cerr << "Error in matrix multiplication - matrices do not produce a valid product!" << std::endl;
		terminaltextcolor(WHITE);
		err::vexit();
	}

	// Set sizes of result matrix C
	const unsigned int C_num_cols = M_num_cols;

	// Declare result matrix
	std::vector<double> C(C_num_cols,0.0);

	// Calculate product
	for(unsigned int j=0;j<C_num_cols;j++){
		for(unsigned int k=0; k<V_num_cols ;k++){
			C[j]+=V[k]*M[k][j];
		}
	}


	// Return
	return C;
}

/// @brief Matrix Transpose Function
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return Matrix B
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///	Reference from Matmul guide at http://www.purplemath.com/modules/mtrxmult3.htm
///=====================================================================================
///
std::vector<std::vector<double> > transpose(
	std::vector<std::vector<double> >& A)
{

	// Get sizes of matrix A
	const unsigned int A_num_rows = A.size();
	const unsigned int A_num_cols = A[0].size();

	// Set sizes of result matrix B
	const unsigned int B_num_rows = A_num_cols;
	const unsigned int B_num_cols = A_num_rows;

	// Declare result matrix
	std::vector<std::vector<double> > B;
	B.resize(B_num_rows);
	for(unsigned int row=0;row<B_num_rows;row++) B[row].resize(B_num_cols,0.0);

	//std::cout << "Transpose is:" << std::endl;
	// C[i][j] = A[j][i]
	// Calculate product
	for(unsigned int i=0;i<B_num_rows;i++){
		//std::cout << "[";
		for(unsigned int j=0;j<B_num_cols;j++){
				B[i][j]=A[j][i];
			//std::cout << B[i][j] << "\t";
		}
		//std::cout << "]"<< std::endl;
	}

	// Return
	return B;
}

/// @brief Matrix Initialisation Function
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    05/02/2011
///
/// @return Matrix C
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///	Reference from Matmul guide at http://www.purplemath.com/modules/mtrxmult3.htm
///=====================================================================================
///
std::vector<std::vector<double> > set_matrix(
	const unsigned int rows,
	const unsigned int cols,
	std::vector<double> & nums)
{

	if(nums.size()!=rows*cols){
		terminaltextcolor(RED);
		std::cerr << "Error in matrix initialisation, incorrect number of elements for matrix!" << std::endl;
		terminaltextcolor(WHITE);
		err::vexit();
	}

	// Declare result matrix
	std::vector<std::vector<double> > C;
	C.resize(rows);
	for(unsigned int row=0;row<rows;row++) C[row].resize(cols,0.0);

	//counter
	unsigned int counter=0;

	//std::cout << "matrix initialised as:" << std::endl;
	for(unsigned int row=0;row<rows;row++){
		//std::cout << "[";
		for(unsigned int col=0;col<cols;col++){
			C[row][col]=nums[counter];
			//std::cout << C[row][col] << "\t";
			counter++;
		}
		//std::cout << "]" << std::endl;
	}

	// Return
	return C;
}

std::vector<std::vector<double> > set_matrix(
	const unsigned int rows,
	const unsigned int cols)
{

	//if(nums.size()!=rows*cols){
	//	std::cerr << "Error in matrix initialisation, incorrect number of elements for matrix!" << std::endl;
	//	err::vexit();
	//}

	// Declare result matrix
	std::vector<std::vector<double> > C;
	C.resize(rows);
	for(unsigned int row=0;row<rows;row++) C[row].resize(cols,0.0);

	//counter
	//unsigned int counter=0;

	//std::cout << "matrix initialised as:" << std::endl;
	for(unsigned int row=0;row<rows;row++){
		//std::cout << "[";
		for(unsigned int col=0;col<cols;col++){
			C[row][col]=0.0;
			//std::cout << C[row][col] << "\t";
			//counter++;
		}
		//std::cout << "]" << std::endl;
	}

	// Return
	return C;
}

void print_matrix( std::vector<std::vector<double> >& A){

	// Get sizes of matrix A
	const unsigned int A_num_rows = A.size();
	const unsigned int A_num_cols = A[0].size();

	for(unsigned int i=0;i<A_num_rows;i++){
		std::cout << "[ ";
		for(unsigned int j=0;j<A_num_cols;j++){
			std::cout << A[i][j] << " ";
		}
		std::cout << "]"<< std::endl;
	}

	// Return
	return;
}

void set_rotational_matrix(
	double ddx,
	double ddy,
	double ddz,
	std::vector< std::vector<double> > & x_rotation_matrix,
	std::vector< std::vector<double> > & y_rotation_matrix,
	std::vector< std::vector<double> > & z_rotation_matrix)
{
	double dx,dy,dz; //change in angle in radians
	double sin_x,sin_y,sin_z,cos_x,cos_y,cos_z;
	double pi=3.14159265358979323846264338327;

	//--------------------------------------------------
	// Initialise varibales
	//--------------------------------------------------

	dx = (ddx/360.0)*2.0*pi;
	dy = (ddy/360.0)*2.0*pi;
	dz = (ddz/360.0)*2.0*pi;

	sin_x = sin(dx);
	cos_x = cos(dx);
	sin_y = sin(dy);
	cos_y = cos(dy);
	sin_z = sin(dz);
	cos_z = cos(dz);

	x_rotation_matrix=vmath::set_matrix(3,3);
	y_rotation_matrix=vmath::set_matrix(3,3);
	z_rotation_matrix=vmath::set_matrix(3,3);

	x_rotation_matrix[0][0] = 1.0;
	x_rotation_matrix[1][0] = 0.0;
	x_rotation_matrix[2][0] = 0.0;
	x_rotation_matrix[0][1] = 0.0;
	x_rotation_matrix[1][1] = cos_x;
	x_rotation_matrix[2][1] = -sin_x;
	x_rotation_matrix[0][2] = 0.0;
	x_rotation_matrix[1][2] = sin_x;
	x_rotation_matrix[2][2] = cos_x;

	y_rotation_matrix[0][0] = cos_y;
	y_rotation_matrix[1][0] = 0.0;
	y_rotation_matrix[2][0] = sin_y;
	y_rotation_matrix[0][1] = 0.0;
	y_rotation_matrix[1][1] = 1.0;
	y_rotation_matrix[2][1] = 0.0;
	y_rotation_matrix[0][2] = -sin_y;
	y_rotation_matrix[1][2] = 0.0;
	y_rotation_matrix[2][2] = cos_y;

	z_rotation_matrix[0][0] = cos_z;
	z_rotation_matrix[1][0] = -sin_z;
	z_rotation_matrix[2][0] = 0.0;
	z_rotation_matrix[0][1] = sin_z;
	z_rotation_matrix[1][1] = cos_z;
	z_rotation_matrix[2][1] = 0.0;
	z_rotation_matrix[0][2] = 0.0;
	z_rotation_matrix[1][2] = 0.0;
	z_rotation_matrix[2][2] = 1.0;

}

// sign functions
double sign(double a){
		if(a<0.0) return -1.0;
		else return 1.0;
}

int sign(int a){
		if(a<0) return -1;
		else return 1;
}

//-----------------------------------------------------
// functions to interpolate values xi,yi -> xj,yj
// with straight line fit
//
//    dx = xj - xi
//    dy = yj - yi
//
//    m = dy/dx
//
//    y = mx + c, solve for c
//
//    c = yi-m*xi
//
double interpolate_m(double xi,double yi,double xj,double yj){

   return (yj-yi)/(xj-xi);

}

double interpolate_c(double xi,double yi,double xj,double yj){

   double m=(yj-yi)/(xj-xi);

   return yi - m*xi;

}

double minimum3(double first, double second, double third){
   double minimum=(std::numeric_limits<double>::max());
   if(first < minimum) minimum = first; // probably always true
   if(second < minimum) minimum = second;
   if(third < minimum) minimum = third;
   return minimum;
}

//--------------------------------------------------------------------------------
// Function to implement linear regression, extracting gradient m and intercept c
//--------------------------------------------------------------------------------
void regression(std::vector<double>& x, std::vector<double>& y, double& m, double& c){

	// check sizes of x and y arrays are the same
	if( x.size() != y.size() ){
	 	std::cerr << "Error in linear regression: x and y array sizes are not the same!" << std::endl;
	 	err::vexit();
	}

	// calculate number of points in data set
	const double num_points = double( x.size() );

	// calculate mean x value
	double sum_x = 0.0;
	for(unsigned int i=0; i < x.size(); i++) sum_x += x[i];

	// calculate mean y value
	double sum_y = 0.0;
	for(unsigned int i=0; i < y.size(); i++) sum_y += y[i];

	// determine sum of squares of x-deviations
	double sum_xx = 0.0;
	for(unsigned int i=0; i< x.size(); i++) sum_xx += x[i]*x[i];

	// determine sum of squares of xy-deviations
	double sum_xy	= 0.0;
	for(unsigned int i=0; i< x.size(); i++) sum_xy += x[i]*y[i];

	// compute gradient
	m = (num_points * sum_xy - sum_x * sum_y) / ( num_points * sum_xx - sum_x * sum_x );

	// compute intercept
	c = ( sum_y - m * sum_x ) / num_points;

	return;

}


} // end of namespcae vmath
