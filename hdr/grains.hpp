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


