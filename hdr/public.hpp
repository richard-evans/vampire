#ifndef PUBLIC_H_
#define PUBLIC_H_
//======================================================================
//                       Global Atomistic Variables
//======================================================================
//#include <vector>
//#include <fstream>
//#include <iostream>
//#include <string>
//#include <cmath>
//#include "mtrand.h"

//#include <valarray>

//using namespace std;






//==========================================================
// Namespace errors
//==========================================================
namespace error_checking
{
	extern bool error_check;
	
}
//==========================================================
// Namespace grains
//==========================================================
namespace grains{

	extern int num_grains;
	
	extern int*  grain_volume_array;
	extern double** grain_coord_array;
	extern double** grain_spin_array;
	extern int* num_assoc_vertices_array;
	extern double** grain_pointx_array;
	extern double** grain_pointy_array;

}



#endif /*PUBLIC_H_*/

