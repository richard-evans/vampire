//====================================================================
//                           set_atom_vars
//   Subroutine to copy newly created system variables to 
//   atom variables for use in integration subroutines
//
//==================================================================== 

#include <iostream>
#include <vector>
#include "atoms.hpp"
#include "material.hpp"
#include "public.hpp"
#include "vmpi.hpp"
#include "create.hpp"

//using namespace atom_variables;
//using namespace material_parameters;
using std::vector;

//==========================================================
// Namespace atom variables
//==========================================================
namespace atoms
{
	//--------------------------
	// Single Variables
	//--------------------------
	int num_atoms;			// Number of atoms in simulation
	int num_neighbours;	   	// Maximum number of neighbours for Hamiltonian/Lattice
	int total_num_neighbours;
	//--------------------------
	// Array Variables
	//--------------------------
	
  vector <double> x_coord_array(0);
	vector <double> y_coord_array(0);
	vector <double> z_coord_array(0);
	vector <int> neighbour_list_array(0);
	vector <int> neighbour_list_start_index(0);
	vector <int> neighbour_list_end_index(0);
	vector <int> type_array(0);
	vector <int> grain_array(0);

	vector <double> x_spin_array(0);
	vector <double> y_spin_array(0);
	vector <double> z_spin_array(0);

	vector <double> x_total_spin_field_array(0);		// Total spin dependent fields
	vector <double> y_total_spin_field_array(0);		// Total spin dependent fields
	vector <double> z_total_spin_field_array(0);		// Total spin dependent fields
	vector <double> x_total_external_field_array(0);	// Total external fields
	vector <double> y_total_external_field_array(0);	// Total external fields
	vector <double> z_total_external_field_array(0);	// Total external fields
	vector <double> x_dipolar_field_array(0);			// Dipolar fields
	vector <double> y_dipolar_field_array(0);			// Dipolar fields
	vector <double> z_dipolar_field_array(0);			// Dipolar fields

}
	
namespace cs{
int set_atom_vars(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <int> > & cneighbourlist){

	// check calling of routine if error checking is activated
	if(error_checking::error_check==true){
		std::cout << "cs::set_atom_vars has been called " << vmpi::my_rank << std::endl;
	}

	//-------------------------------------------------
	// Set number of atoms
	//-------------------------------------------------

	atoms::num_atoms = catom_array.size();
	std::cout << "rank:\t" << vmpi::my_rank << "\tnum atoms:\t" << atoms::num_atoms-vmpi::num_halo_atoms << std::endl; 

	atoms::x_coord_array.resize(atoms::num_atoms,0);
	atoms::y_coord_array.resize(atoms::num_atoms,0);
	atoms::z_coord_array.resize(atoms::num_atoms,0);

	atoms::x_spin_array.resize(atoms::num_atoms,0.0);
	atoms::y_spin_array.resize(atoms::num_atoms,0.0);
	atoms::z_spin_array.resize(atoms::num_atoms,1.0);

	atoms::type_array.resize(atoms::num_atoms,0);
	
	atoms::x_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::y_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::z_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::x_total_external_field_array.resize(atoms::num_atoms,0.0);	
	atoms::y_total_external_field_array.resize(atoms::num_atoms,0.0);	
	atoms::z_total_external_field_array.resize(atoms::num_atoms,0.0);	

	for(int atom=0;atom<atoms::num_atoms;atom++){
		
		atoms::x_coord_array[atom] = catom_array[atom].x;
		atoms::y_coord_array[atom] = catom_array[atom].y;
		atoms::z_coord_array[atom] = catom_array[atom].z;
		
		atoms::type_array[atom] = catom_array[atom].material;
	}

	//===========================================================
	// Create 1-D neighbourlist
	//===========================================================

	//-------------------------------------------------
	//	Calculate total number of neighbours
	//-------------------------------------------------
	int counter = 0;
	
	for(int atom=0;atom<atoms::num_atoms;atom++){
		counter+=cneighbourlist[atom].size();
	}
	
	atoms::total_num_neighbours = counter;
	
	atoms::neighbour_list_array.resize(atoms::total_num_neighbours,0);
	atoms::neighbour_list_start_index.resize(atoms::num_atoms,0);
	atoms::neighbour_list_end_index.resize(atoms::num_atoms,0);

	//	Populate 1D neighbourlist and index arrays
	counter = 0;
	for(int atom=0;atom<atoms::num_atoms;atom++){
		//std::cout << atom << ": ";
		// Set start index
		atoms::neighbour_list_start_index[atom]=counter;
		for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){
			atoms::neighbour_list_array[counter] = cneighbourlist[atom][nn];
			//std::cout << cneighbourlist[atom][nn] << " ";
			counter++;
		}
		//std::cout << std::endl;
		// Set end index
		atoms::neighbour_list_end_index[atom]=counter-1;
	}
	
	return EXIT_SUCCESS;
}

} // End of cs namespace
