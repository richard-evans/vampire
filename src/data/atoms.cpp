#include "atoms.hpp"

#include <vector>

//==========================================================
// Namespace atom variables
//==========================================================
namespace atoms{
	//--------------------------
	// Single Variables
	//--------------------------
	int num_atoms;			// Number of atoms in simulation
	int num_neighbours;	   	// Maximum number of neighbours for Hamiltonian/Lattice
	int total_num_neighbours;
	int exchange_type;
	//--------------------------
	// Array Variables
	//--------------------------
	
	std::vector <double> x_coord_array(0);
	std::vector <double> y_coord_array(0);
	std::vector <double> z_coord_array(0);
	std::vector <int> neighbour_list_array(0);
	std::vector <int> neighbour_interaction_type_array(0);
	std::vector <int> neighbour_list_start_index(0);
	std::vector <int> neighbour_list_end_index(0);
	std::vector <int> type_array(0);
	std::vector <int> category_array(0);
	std::vector <int> grain_array(0);
	std::vector <int> cell_array(0);

	std::vector <double> x_spin_array(0);
	std::vector <double> y_spin_array(0);
	std::vector <double> z_spin_array(0);

	std::vector <double> x_total_spin_field_array(0);		// Total spin dependent fields
	std::vector <double> y_total_spin_field_array(0);		// Total spin dependent fields
	std::vector <double> z_total_spin_field_array(0);		// Total spin dependent fields
	std::vector <double> x_total_external_field_array(0);	// Total external fields
	std::vector <double> y_total_external_field_array(0);	// Total external fields
	std::vector <double> z_total_external_field_array(0);	// Total external fields
	std::vector <double> x_dipolar_field_array(0);			// Dipolar fields
	std::vector <double> y_dipolar_field_array(0);			// Dipolar fields
	std::vector <double> z_dipolar_field_array(0);			// Dipolar fields
	
	std::vector <zval_t> i_exchange_list(0);
	std::vector <zvec_t> v_exchange_list(0);
	std::vector <zten_t> t_exchange_list(0);

}
