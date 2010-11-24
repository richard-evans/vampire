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
	extern int size;
	extern int update_rate;
	extern int update_counter;

	extern std::vector <int> num_atoms_in_cell;

	extern std::vector <double> x_coord_array;
	extern std::vector <double> y_coord_array;
	extern std::vector <double> z_coord_array;

	extern std::vector <double> x_mag_array;
	extern std::vector <double> y_mag_array;
	extern std::vector <double> z_mag_array;
	
	extern std::vector <double> x_field_array;
	extern std::vector <double> y_field_array;
	extern std::vector <double> z_field_array;

	extern int initialise();
	extern int mag();
	extern int output_mag(std::ofstream&);
}

#endif /*CELLS_H_*/


