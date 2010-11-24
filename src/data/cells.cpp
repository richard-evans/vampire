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

// calculate the magnetisation in each cell
int cells::mag() {
  using namespace atoms;

  for(int i=0;i<num_atoms;++i) {
    // lookup which cell atom is in
    int cell = cell_array[i];
    int type = type_array[i];
    const double mus = mp::material[type].mu_s_SI;

    //
    x_mag_array[cell] += x_spin_array[i]*mus;
    y_mag_array[cell] += y_spin_array[i]*mus;
    z_mag_array[cell] += z_spin_array[i]*mus;
  }

  return EXIT_SUCCESS;
}
