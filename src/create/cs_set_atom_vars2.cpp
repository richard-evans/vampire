//====================================================================
//                           set_atom_vars
//   Subroutine to copy newly created system variables to 
//   atom variables for use in integration subroutines
//
//==================================================================== 

#include <iostream>
#include <vector>

#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "vmpi.hpp"


//using namespace atom_variables;
//using namespace material_parameters;
	
namespace cs{
int set_atom_vars(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <neighbour_t> > & cneighbourlist){

	// check calling of routine if error checking is activated
	if(err::check==true){
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
	atoms::category_array.resize(atoms::num_atoms,0);
	atoms::grain_array.resize(atoms::num_atoms,0);
	atoms::cell_array.resize(atoms::num_atoms,0);
	
	atoms::x_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::y_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::z_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::x_total_external_field_array.resize(atoms::num_atoms,0.0);	
	atoms::y_total_external_field_array.resize(atoms::num_atoms,0.0);	
	atoms::z_total_external_field_array.resize(atoms::num_atoms,0.0);	
	atoms::x_dipolar_field_array.resize(atoms::num_atoms,0.0);	
	atoms::y_dipolar_field_array.resize(atoms::num_atoms,0.0);	
	atoms::z_dipolar_field_array.resize(atoms::num_atoms,0.0);	

	for(int atom=0;atom<atoms::num_atoms;atom++){
		
		atoms::x_coord_array[atom] = catom_array[atom].x;
		atoms::y_coord_array[atom] = catom_array[atom].y;
		atoms::z_coord_array[atom] = catom_array[atom].z;
		
		atoms::type_array[atom] = catom_array[atom].material;
		atoms::category_array[atom] = catom_array[atom].lh_category;
		//std::cout << atom << " grain: " << catom_array[atom].grain << std::endl;
		atoms::grain_array[atom] = catom_array[atom].grain;

		// initialise atomic spin positions
		int mat=atoms::type_array[atom];
		double sx,sy,sz; // spins 
		if(mp::material[mat].random_spins==true){
			sx=2.0*mtrandom::grnd()-1.0;
			sy=2.0*mtrandom::grnd()-1.0;
			sz=2.0*mtrandom::grnd()-1.0;
		}
		else{
			sx=mp::material[mat].initial_spin[0];
			sy=mp::material[mat].initial_spin[1];
			sz=mp::material[mat].initial_spin[2];
		}
		// now normalise spins
		double modS=1.0/sqrt(sx*sx + sy*sy + sz*sz);
		atoms::x_spin_array[atom]=sx*modS;
		atoms::y_spin_array[atom]=sy*modS;
		atoms::z_spin_array[atom]=sz*modS;
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
	atoms::neighbour_interaction_type_array.resize(atoms::total_num_neighbours,0);
	atoms::neighbour_list_start_index.resize(atoms::num_atoms,0);
	atoms::neighbour_list_end_index.resize(atoms::num_atoms,0);

	//	Populate 1D neighbourlist and index arrays
	counter = 0;
	for(int atom=0;atom<atoms::num_atoms;atom++){
		//std::cout << atom << ": ";
		// Set start index
		atoms::neighbour_list_start_index[atom]=counter;
		for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){
			atoms::neighbour_list_array[counter] = cneighbourlist[atom][nn].nn;
			atoms::neighbour_interaction_type_array[counter] = cneighbourlist[atom][nn].i;
			//std::cout << cneighbourlist[atom][nn] << " ";
			counter++;
		}
		//std::cout << std::endl;
		// Set end index
		atoms::neighbour_list_end_index[atom]=counter-1;
	}
	
	// condense interaction list
	atoms::exchange_type=unit_cell.exchange_type;
	
	// temporary class variables
	zval_t tmp_zval;
	zvec_t tmp_zvec;
	zten_t tmp_zten;
	
	switch(atoms::exchange_type){
		case -1:
			// unroll material calculations
			std::cout << "Unrolled exchange template requires " << 1.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			atoms::i_exchange_list.reserve(atoms::neighbour_list_array.size());
			// loop over all interactions
			for(int atom=0;atom<atoms::num_atoms;atom++){
				const int imaterial=atoms::type_array[atom];
				for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int jmaterial=atoms::type_array[natom];
					atoms::i_exchange_list.push_back(tmp_zval);
					atoms::i_exchange_list[atom].Jij= mp::material[imaterial].Jij_matrix[jmaterial];
				}
			}
			// now set exchange type to normal isotropic case
			atoms::exchange_type=0;
			break;
		case 0:
			std::cout << "Unrolled exchange template requires " << 1.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::i_exchange_list.reserve(unit_cell.interaction.size());
			for(int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::i_exchange_list.push_back(tmp_zval);
				atoms::i_exchange_list[i].Jij=unit_cell.interaction[i].Jij[0][0]*1.0e-21/mp::material[imat].mu_s_SI;
			}
			break;
		case 1:
			std::cout << "Unrolled exchange template requires " << 3.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::v_exchange_list.reserve(unit_cell.interaction.size());
			for(int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::v_exchange_list.push_back(tmp_zvec);
				atoms::v_exchange_list[i].Jij[0]=unit_cell.interaction[i].Jij[0][0]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::v_exchange_list[i].Jij[1]=unit_cell.interaction[i].Jij[1][1]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::v_exchange_list[i].Jij[2]=unit_cell.interaction[i].Jij[2][2]*1.0e-21/mp::material[imat].mu_s_SI;
			}
			break;
		case 2:
			std::cout << "Unrolled exchange template requires " << 9.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::t_exchange_list.reserve(unit_cell.interaction.size());
			for(int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::t_exchange_list.push_back(tmp_zten);
				
				atoms::t_exchange_list[i].Jij[0][0]=unit_cell.interaction[i].Jij[0][0]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[0][1]=unit_cell.interaction[i].Jij[0][1]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[0][2]=unit_cell.interaction[i].Jij[0][2]*1.0e-21/mp::material[imat].mu_s_SI;

				atoms::t_exchange_list[i].Jij[1][0]=unit_cell.interaction[i].Jij[1][0]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[1][1]=unit_cell.interaction[i].Jij[1][1]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[1][2]=unit_cell.interaction[i].Jij[1][2]*1.0e-21/mp::material[imat].mu_s_SI;

				atoms::t_exchange_list[i].Jij[2][0]=unit_cell.interaction[i].Jij[2][0]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[2][1]=unit_cell.interaction[i].Jij[2][1]*1.0e-21/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[2][2]=unit_cell.interaction[i].Jij[2][2]*1.0e-21/mp::material[imat].mu_s_SI;
			}
			break;
		default:
			std::cerr << "Error! - Unknown unit cell exchange type " << atoms::exchange_type << "; unable to unroll exchenge template. Exiting" << std::endl;
			err::vexit();
			break;
	}
	
	// now remove unit cell interactions data
	unit_cell.interaction.resize(0);
	unit_cell.atom.resize(0);
	

	return EXIT_SUCCESS;
}

} // End of cs namespace
