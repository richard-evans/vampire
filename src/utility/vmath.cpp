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

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace vmath{

bool point_in_polygon(double x, double y,double *polyX, double *polyY,int polySides) {
	//========================================================================================================
	//		 						Function to decide if point is within polygon
	//
	//														Version 1.0
	//
	//												R F Evans 15/07/2009
	//========================================================================================================

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
		std::cerr << "Error in matrix multiplication - matrices do not produce a valid product!" << std::endl;
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
		std::cerr << "Error in matrix initialisation, incorrect number of elements for matrix!" << std::endl;
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
		  
// sign functions
double sign(double a){
		if(a<0.0) return -1.0;
		else return 1.0;
}

int sign(int a){
		if(a<0) return -1;
		else return 1;
}

} // end of namespcae vmath

