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

int cells::mag() {
  using namespace atoms;

  // calulate magnetisation in each cell
  for(int i=0;i<num_atoms;++i) {
    int cell = cell_array[i];
    int type = type_array[i];
    const double mus = mp::material[type].mu_s_SI;

    //
    x_mag_array[cell] += x_spin_array[i]*mus;
    y_mag_array[cell] += y_spin_array[i]*mus;
    z_mag_array[cell] += z_spin_array[i]*mus;
  }

#ifdef MPICF
  // reduce (sum) arrays to root node
  MPI::COMM_WORLD.Reduce( &(x_mag_array[0]), &(x_mag_array[0]), num_cells,
      MPI_DOUBLE, MPI_SUM, 0 );
  MPI::COMM_WORLD.Reduce( &(y_mag_array[0]), &(y_mag_array[0]), num_cells,
      MPI_DOUBLE, MPI_SUM, 0 );
  MPI::COMM_WORLD.Reduce( &(z_mag_array[0]), &(z_mag_array[0]), num_cells,
      MPI_DOUBLE, MPI_SUM, 0 );

  // broadcast result of reduction from root node
  MPI::COMM_WORLD.Bcast( &(x_mag_array[0]), num_cells, MPI_DOUBLE, 0 );
  MPI::COMM_WORLD.Bcast( &(y_mag_array[0]), num_cells, MPI_DOUBLE, 0 );
  MPI::COMM_WORLD.Bcast( &(z_mag_array[0]), num_cells, MPI_DOUBLE, 0 );
#endif

  return EXIT_SUCCESS;
}
