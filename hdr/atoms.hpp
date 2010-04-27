#ifndef ATOMS_H_
#define ATOMS_H_

#include <string>
#include <vector>

//======================================================================
//                       Global Atomistic Variables
//======================================================================
namespace atoms
{
	//--------------------------
	// Single Variables
	//--------------------------
	extern int num_atoms;			// Number of atoms in simulation
	extern int num_neighbours;	   	// Maximum number of neighbours for Hamiltonian/Lattice
	extern int total_num_neighbours;// Total number of neighbours for system
	//--------------------------
	// Array Variables
	//--------------------------

	extern std::vector <double> x_coord_array;
	extern std::vector <double> y_coord_array;
	extern std::vector <double> z_coord_array;
	extern std::vector <int> neighbour_list_array;
	extern std::vector <int> neighbour_list_start_index;
	extern std::vector <int> neighbour_list_end_index;
	extern std::vector <int> type_array;
	extern std::vector <int> grain_array;

	extern std::vector <double> x_spin_array;
	extern std::vector <double> y_spin_array;
	extern std::vector <double> z_spin_array;

	extern std::vector <double> x_total_spin_field_array;		// Total spin dependent fields
	extern std::vector <double> y_total_spin_field_array;		// Total spin dependent fields
	extern std::vector <double> z_total_spin_field_array;		// Total spin dependent fields
	extern std::vector <double> x_total_external_field_array;	// Total external fields
	extern std::vector <double> y_total_external_field_array;	// Total external fields
	extern std::vector <double> z_total_external_field_array;	// Total external fields
	extern std::vector <double> x_dipolar_field_array;			// Dipolar fields
	extern std::vector <double> y_dipolar_field_array;			// Dipolar fields
	extern std::vector <double> z_dipolar_field_array;			// Dipolar fields
}


#endif /*ATOMS_H_*/
