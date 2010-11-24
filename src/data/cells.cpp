#include "atoms.hpp"  
#include "cells.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include <cmath>
#include <iostream>

namespace cells{
	
	int num_cells;
	int size;
	int update_rate;
	int update_counter;

	std::vector <int> num_atoms_in_cell;

	std::vector <double> x_coord_array;
	std::vector <double> y_coord_array;
	std::vector <double> z_coord_array;

	std::vector <double> x_mag_array;
	std::vector <double> y_mag_array;
	std::vector <double> z_mag_array;
	
	std::vector <double> x_field_array;
	std::vector <double> y_field_array;
	std::vector <double> z_field_array;
	
	int initialise();
	int mag();
	int output_mag(std::ofstream&);
	
} // End of namespace cells
