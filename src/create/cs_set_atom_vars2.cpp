//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//====================================================================
//                           set_atom_vars
//   Subroutine to copy newly created system variables to
//   atom variables for use in integration subroutines
//
//====================================================================

#include <iostream>
#include <vector>

#include "dipole.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "grains.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"


//using namespace atom_variables;
//using namespace material_parameters;

namespace cs{
int set_atom_vars(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <neighbour_t> > & cneighbourlist){

	// check calling of routine if error checking is activated
	if(err::check==true){
		std::cout << "cs::set_atom_vars has been called " << vmpi::my_rank << std::endl;
	}

	// Set number of atoms
	atoms::num_atoms = catom_array.size(); // core and boundary spins in mpi mode

	zlog << zTs() << "Number of atoms generated on rank " << vmpi::my_rank << ": " << atoms::num_atoms-vmpi::num_halo_atoms << std::endl;
	zlog << zTs() << "Memory required for copying to performance array on rank " << vmpi::my_rank << ": " << 19.0*double(atoms::num_atoms)*8.0/1.0e6 << " MB RAM"<< std::endl;

   // Save number of non-magnetic atoms
   atoms::num_non_magnetic_atoms = cs::non_magnetic_atoms_array.size();

   atoms::x_coord_array.resize(atoms::num_atoms,0.0);
   atoms::y_coord_array.resize(atoms::num_atoms,0.0);
   atoms::z_coord_array.resize(atoms::num_atoms,0.0);

	atoms::x_spin_array.resize(atoms::num_atoms,0.0);
	atoms::y_spin_array.resize(atoms::num_atoms,0.0);
	atoms::z_spin_array.resize(atoms::num_atoms,1.0);
   atoms::m_spin_array.resize(atoms::num_atoms,0.0);

   atoms::type_array.resize(     atoms::num_atoms,0);
   atoms::category_array.resize( atoms::num_atoms,0);
   atoms::grain_array.resize(    atoms::num_atoms,0);
   atoms::cell_array.resize(     atoms::num_atoms,0);

	atoms::x_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::y_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::z_total_spin_field_array.resize(atoms::num_atoms,0.0);
	atoms::x_total_external_field_array.resize(atoms::num_atoms,0.0);
	atoms::y_total_external_field_array.resize(atoms::num_atoms,0.0);
	atoms::z_total_external_field_array.resize(atoms::num_atoms,0.0);
//	atoms::x_dipolar_field_array.resize(atoms::num_atoms,0.0);
//	atoms::y_dipolar_field_array.resize(atoms::num_atoms,0.0);
//	atoms::z_dipolar_field_array.resize(atoms::num_atoms,0.0);
   // Resize to zero atoms_dipolar_field-x,y,z
	dipole::atom_dipolar_field_array_x.resize(atoms::num_atoms,0.0);
	dipole::atom_dipolar_field_array_y.resize(atoms::num_atoms,0.0);
	dipole::atom_dipolar_field_array_z.resize(atoms::num_atoms,0.0);
   // Resize to zero atoms_mu0demag_field-x,y,z
	dipole::atom_mu0demag_field_array_x.resize(atoms::num_atoms,0.0);
	dipole::atom_mu0demag_field_array_y.resize(atoms::num_atoms,0.0);
	dipole::atom_mu0demag_field_array_z.resize(atoms::num_atoms,0.0);

   // Set custom RNG for spin initialisation
   MTRand random_spin_rng;
   random_spin_rng.seed(123456+vmpi::my_rank);

	for(int atom=0;atom<atoms::num_atoms;atom++){

		atoms::x_coord_array[atom] = catom_array[atom].x;
		atoms::y_coord_array[atom] = catom_array[atom].y;
		atoms::z_coord_array[atom] = catom_array[atom].z;

		atoms::type_array[atom] = catom_array[atom].material;
		atoms::category_array[atom] = catom_array[atom].lh_category;
		//std::cout << atom << " grain: " << catom_array[atom].grain << std::endl;
		atoms::grain_array[atom] = catom_array[atom].grain;

		// initialise atomic spin positions
      // Use a normalised gaussian for uniform distribution on a unit sphere
		int mat=atoms::type_array[atom];
		double sx,sy,sz; // spins
		if(mp::material[mat].random_spins==true){
         sx=mtrandom::gaussianc(random_spin_rng);
         sy=mtrandom::gaussianc(random_spin_rng);
         sz=mtrandom::gaussianc(random_spin_rng);
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
      atoms::m_spin_array[atom]=mp::material[mat].mu_s_SI/9.27400915e-24;
	}

	//===========================================================
	// Create 1-D neighbourlist
	//===========================================================

	zlog << zTs() << "Memory required for creation of 1D neighbour list on rank " << vmpi::my_rank << ": ";
	zlog << (2.0*double(atoms::num_atoms)+2.0*double(atoms::total_num_neighbours))*8.0/1.0e6 << " MB RAM"<< std::endl;

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
			if(cneighbourlist[atom][nn].nn > atoms::num_atoms){
				terminaltextcolor(RED);
				std::cerr << "Fatal Error - neighbour " << cneighbourlist[atom][nn].nn <<" is out of valid range 0-"
				<< atoms::num_atoms << " on rank " << vmpi::my_rank << std::endl;
				std::cerr << "Atom " << atom << " of MPI type " << catom_array[atom].mpi_type << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}

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
			std::cout << "Using generic form of exchange interaction with " << unit_cell.interaction.size() << " total interactions." << std::endl;
			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			atoms::i_exchange_list.reserve(atoms::neighbour_list_array.size());
			// loop over all interactions
			for(int atom=0;atom<atoms::num_atoms;atom++){
				const int imaterial=atoms::type_array[atom];
				for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
					const int natom = atoms::neighbour_list_array[nn];
					const int jmaterial=atoms::type_array[natom];
					atoms::i_exchange_list.push_back(tmp_zval);
               // get unit cell interaction id
               int i = atoms::neighbour_interaction_type_array[nn];
               atoms::i_exchange_list[nn].Jij=unit_cell.interaction[i].Jij[0][0]*mp::material[imaterial].Jij_matrix[jmaterial][0];
					// reset interation id to neighbour number - causes segfault if nn out of range
					atoms::neighbour_interaction_type_array[nn]=nn;
				}
			}
			// now set exchange type to normal isotropic case
			atoms::exchange_type=0;
			break;
		case 0:
			std::cout << "Using isotropic form of exchange interaction with " << unit_cell.interaction.size() << " total interactions." << std::endl;
			zlog << zTs() << "Unrolled exchange template requires " << 1.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::i_exchange_list.reserve(unit_cell.interaction.size());
			for(unsigned int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::i_exchange_list.push_back(tmp_zval);
				atoms::i_exchange_list[i].Jij=-unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
			}
			break;
		case 1:
			std::cout << "Using vectorial form of exchange interaction with " << unit_cell.interaction.size() << " total interactions." << std::endl;
			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::v_exchange_list.reserve(unit_cell.interaction.size());
			for(unsigned int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::v_exchange_list.push_back(tmp_zvec);
				atoms::v_exchange_list[i].Jij[0]=-unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
				atoms::v_exchange_list[i].Jij[1]=-unit_cell.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
				atoms::v_exchange_list[i].Jij[2]=-unit_cell.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
			}
			break;
		case 2:
			std::cout << "Using tensorial form of exchange interaction with " << unit_cell.interaction.size() << " total interactions." << std::endl;
			zlog << zTs() << "Unrolled exchange template requires " << 9.0*double(unit_cell.interaction.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
			// unroll isotopic interactions
			atoms::t_exchange_list.reserve(unit_cell.interaction.size());
			for(unsigned int i=0;i<unit_cell.interaction.size();i++){
				int iatom = unit_cell.interaction[i].i;
				int imat = unit_cell.atom[iatom].mat;
				atoms::t_exchange_list.push_back(tmp_zten);

				atoms::t_exchange_list[i].Jij[0][0]=-unit_cell.interaction[i].Jij[0][0]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[0][1]=-unit_cell.interaction[i].Jij[0][1]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[0][2]=-unit_cell.interaction[i].Jij[0][2]/mp::material[imat].mu_s_SI;

				atoms::t_exchange_list[i].Jij[1][0]=-unit_cell.interaction[i].Jij[1][0]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[1][1]=-unit_cell.interaction[i].Jij[1][1]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[1][2]=-unit_cell.interaction[i].Jij[1][2]/mp::material[imat].mu_s_SI;

				atoms::t_exchange_list[i].Jij[2][0]=-unit_cell.interaction[i].Jij[2][0]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[2][1]=-unit_cell.interaction[i].Jij[2][1]/mp::material[imat].mu_s_SI;
				atoms::t_exchange_list[i].Jij[2][2]=-unit_cell.interaction[i].Jij[2][2]/mp::material[imat].mu_s_SI;
			}
			break;

      case 3: // normalised vectorial exchange
   			// unroll material calculations
   			std::cout << "Using vectorial form of exchange interaction with " << unit_cell.interaction.size() << " total interactions." << std::endl;
   			zlog << zTs() << "Unrolled exchange template requires " << 3.0*double(atoms::neighbour_list_array.size())*double(sizeof(double))*1.0e-6 << "MB RAM" << std::endl;
   			atoms::v_exchange_list.reserve(atoms::neighbour_list_array.size());
   			// loop over all interactions
   			for(int atom=0;atom<atoms::num_atoms;atom++){
   				const int imaterial=atoms::type_array[atom];
   				for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
   					const int natom = atoms::neighbour_list_array[nn];
   					const int jmaterial=atoms::type_array[natom];
   					atoms::v_exchange_list.push_back(tmp_zvec);
                  // get unit cell interaction id
                  int i = atoms::neighbour_interaction_type_array[nn];
                  atoms::v_exchange_list[nn].Jij[0]=unit_cell.interaction[i].Jij[0][0]*mp::material[imaterial].Jij_matrix[jmaterial][0];
                  atoms::v_exchange_list[nn].Jij[1]=unit_cell.interaction[i].Jij[1][1]*mp::material[imaterial].Jij_matrix[jmaterial][1];
                  atoms::v_exchange_list[nn].Jij[2]=unit_cell.interaction[i].Jij[2][2]*mp::material[imaterial].Jij_matrix[jmaterial][2];
   					// reset interation id to neighbour number - causes segfault if nn out of range
   					atoms::neighbour_interaction_type_array[nn]=nn;
   				}
   			}
   			// now set exchange type to normal vectorial case
   			atoms::exchange_type=1;
   			break;
		default:
			terminaltextcolor(RED);
			std::cerr << "Error! - Unknown unit cell exchange type " << atoms::exchange_type << "; unable to unroll exchenge template. Exiting" << std::endl;
			terminaltextcolor(WHITE);
			err::vexit();
			break;
	}



   // now remove unit cell interactions data
   unit_cell.interaction.resize(0);

   // Save number of atoms in unit cell first
   cells::num_atoms_in_unit_cell=unit_cell.atom.size();
   //std::cout << "\t\t" << unit_cell.atom.size() << "\t" << cells::num_atoms_in_unit_cell << std::endl;
   unit_cell.atom.resize(0);

   // Now nuke generation vectors to free memory NOW
   std::vector<cs::catom_t> zerov;
   std::vector<std::vector <neighbour_t> > zerovv;
   catom_array.swap(zerov);
   cneighbourlist.swap(zerovv);



   return EXIT_SUCCESS;

}

} // End of cs namespace
