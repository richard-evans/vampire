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
}