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

	//-------------------------------------------------
	// Set number of atoms
	//-------------------------------------------------

	atoms::num_atoms = catom_array.size();
	zlog << zTs() << "Number of atoms generated on rank " << vmpi::my_rank << ": " << atoms::num_atoms-vmpi::num_halo_atoms << std::endl;
	zlog << zTs() << "Memory required for copying to performance array on rank " << vmpi::my_rank << ": " << 19.0*double(atoms::num_atoms)*8.0/1.0e6 << " MB RAM"<< std::endl;

	atoms::x_coord_array.resize(atoms::num_atoms,0);
	atoms::y_coord_array.resize(atoms::num_atoms,0);
	atoms::z_coord_array.resize(atoms::num_atoms,0);

	atoms::x_spin_array.resize(atoms::num_atoms,0.0);
	atoms::y_spin_array.resize(atoms::num_atoms,0.0);
	atoms::z_spin_array.resize(atoms::num_atoms,1.0);
   atoms::m_spin_array.resize(atoms::num_atoms,0.0);

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

	// initialise surface threshold if not overidden by input file
	if(sim::surface_anisotropy_threshold==123456789) sim::surface_anisotropy_threshold=unit_cell.surface_threshold;

	//-------------------------------------------------
	//	Optionally set up surface anisotropy
	//-------------------------------------------------

	// create temporary array for storing surface threshold
	std::vector<unsigned int> surface_anisotropy_threshold_array(atoms::num_atoms, sim::surface_anisotropy_threshold);
	// if using native (local) surface threshold then repopulate threshold array
	if(sim::NativeSurfaceAnisotropyThreshold){
	  zlog << zTs() << "Identifying surface atoms using native (site dependent) threshold." << std::endl;
	  for(int atom=0;atom<atoms::num_atoms;atom++){
			unsigned int atom_uc_id=catom_array.at(atom).uc_id;
			surface_anisotropy_threshold_array.at(atom)=unit_cell.atom.at(atom_uc_id).ni;
		}
	}
	else zlog << zTs() << "Identifying surface atoms using global threshold value of " << sim::surface_anisotropy_threshold << "." << std::endl;

   //--------------------------------------------------------------------------------------------
   // Determine nearest neighbour interactions from unit cell data for a single unit cell
   //--------------------------------------------------------------------------------------------

   // vector to store interactions withn range
   std::vector<bool> nn_interaction(cs::unit_cell.interaction.size(),false);

   // save nn_distance for performance
   const double rsq=sim::nearest_neighbour_distance*sim::nearest_neighbour_distance;

   // Get unit cell size
   const double ucdx=cs::unit_cell.dimensions[0];
   const double ucdy=cs::unit_cell.dimensions[1];
   const double ucdz=cs::unit_cell.dimensions[2];

   // loop over all interactions in unit cell
   for(unsigned int itr=0;itr<cs::unit_cell.interaction.size();itr++){

      // get distance to neighbouring unit cell in unit cells
      double nndx=double(cs::unit_cell.interaction[itr].dx);
      double nndy=double(cs::unit_cell.interaction[itr].dy);
      double nndz=double(cs::unit_cell.interaction[itr].dz);

      // load positions of i and j atoms to temporary coordinates and convert to angstroms
      double ix=(unit_cell.atom[cs::unit_cell.interaction[itr].i].x)*ucdx;
      double iy=(unit_cell.atom[cs::unit_cell.interaction[itr].i].y)*ucdy;
      double iz=(unit_cell.atom[cs::unit_cell.interaction[itr].i].z)*ucdz;
      double jx=(unit_cell.atom[cs::unit_cell.interaction[itr].j].x + nndx)*ucdx;
      double jy=(unit_cell.atom[cs::unit_cell.interaction[itr].j].y + nndy)*ucdy;
      double jz=(unit_cell.atom[cs::unit_cell.interaction[itr].j].z + nndz)*ucdz;

      // calculate reduced coordinates
      double dx=jx-ix;
      double dy=jy-iy;
      double dz=jz-iz;

      // calculate interaction range and check if less than nn distance
      const double range = (dx*dx + dy*dy + dz*dz);
      if(range <=rsq) nn_interaction[itr]=true;
   }

   //------------------------------------------------------------
   // Identify all nearest neighbour interactions in system
   //
   // Nearest neighbour list is a subset of full neighbour list,
   // and so everything is derived from that.
   //------------------------------------------------------------

   // vector to identify all nearest neighbour interactions
   std::vector <std::vector <bool> > nearest_neighbour_interactions_list(atoms::num_atoms);

   // loop over all atoms
   for(int atom=0;atom<atoms::num_atoms;atom++){

      // set all interactions for atom as non-nearest neighbour by default
      nearest_neighbour_interactions_list[atom].resize(cneighbourlist[atom].size(),false);

      // loop over all interactions for atom
      for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){

         // get interaction type (same as unit cell interaction id)
         unsigned int id = cneighbourlist[atom][nn].i;

         // Ensure valid interaction id
         if(id>nn_interaction.size()){
            std::cout << "Error: invalid interaction id " << id << " is greater than number of interactions in unit cell " << nn_interaction.size() << ". Exiting" << std::endl;
            zlog << zTs() << "Error: invalid interaction id " << id << " is greater than number of interactions in unit cell " << nn_interaction.size() << ". Exiting" << std::endl;
            err::vexit();
         }

         // set mask to true or false for non-fully coordinated atoms in the bulk
         nearest_neighbour_interactions_list[atom][nn]=nn_interaction.at(id);
      }
   }

   //----------------------------------------------------------------------------------------
   // Identify atoms with less than full nearest neighbour coordination
   //----------------------------------------------------------------------------------------

   // Track total number of surface atoms and total nearest neighbour interactions
   unsigned int num_surface_atoms=0;
   unsigned int total_num_surface_nn=0;

   // Resize surface atoms mask and initialise to false
   atoms::surface_array.resize(atoms::num_atoms, false);

   // Loop over all *local* atoms
   for(int atom=0;atom<atoms::num_atoms;atom++){

      // Check for local MPI atoms only
      if(catom_array[atom].mpi_type!=2){

         // Initialise counter for number of nearest neighbour interactions
         unsigned int nnn_int=0;

         // Loop over all interactions to determine number of nearest neighbour interactions
         for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){

            // If interaction is nn, increment counter
            if(nearest_neighbour_interactions_list[atom][nn]) nnn_int++;

         }

         // check for atoms with < threshold number of nearest neighbours
         if(nnn_int<surface_anisotropy_threshold_array.at(atom)){
            atoms::surface_array[atom]=true;
            num_surface_atoms++;
            total_num_surface_nn+=nnn_int;
         }
      }
   }

   // Output statistics to log file
   zlog << zTs() << num_surface_atoms << " surface atoms found." << std::endl;

   //----------------------------------------------------------------
   // If surface anisotropy is enabled, calculate necessary data
   //----------------------------------------------------------------
   if(sim::surface_anisotropy==true){

      zlog << zTs() << "Using surface anisotropy for atoms with < threshold number of neighbours." << std::endl;

      atoms::surface_array.resize(atoms::num_atoms);
      atoms::nearest_neighbour_list_si.resize(atoms::num_atoms);
      atoms::nearest_neighbour_list_ei.resize(atoms::num_atoms);
      atoms::nearest_neighbour_list.reserve(total_num_surface_nn);

      // counter for index arrays
      int counter=0;

      //	Populate surface atom and 1D nearest neighbour list and index arrays
      for(int atom=0;atom<atoms::num_atoms;atom++){

         // initialise 1D index arrays
         atoms::nearest_neighbour_list_si[atom]=0;
         atoms::nearest_neighbour_list_ei[atom]=0;

         //std::cout << "Atom " << atom << std::endl;
         //std::cout << "--------------------------------------------" << std::endl;

         // Only calculate parameters for atoms with less than full nn coordination
         if(atoms::surface_array[atom]){

            // Set start index
            atoms::nearest_neighbour_list_si[atom]=counter;

            // loop over all neighbours
            for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){

               // only add nearest neighbours to list
               if(nearest_neighbour_interactions_list[atom][nn]==true){

                  // add interaction to 1D list
                  atoms::nearest_neighbour_list.push_back(cneighbourlist[atom][nn].nn);

                  // get atomic position vector i->j
                  double eij[3]={cneighbourlist[atom][nn].vx,cneighbourlist[atom][nn].vy,cneighbourlist[atom][nn].vz};

                  // normalise to unit vector
                  double invrij=1.0/sqrt(eij[0]*eij[0]+eij[1]*eij[1]+eij[2]*eij[2]);

                  atoms::eijx.push_back(eij[0]*invrij);
                  atoms::eijy.push_back(eij[1]*invrij);
                  atoms::eijz.push_back(eij[2]*invrij);

                  //int natom = cneighbourlist[atom][nn].nn;
                  //std::cout << "nn_id: " << nn << " j: " << cneighbourlist[atom][nn].nn << "\trange: " << 1.0/invrij << " ";
                  //std::cout << "eij: " << eij[0] << " " << eij[1] << " " << eij[2] << "\tatomi: ";
                  //std::cout << catom_array[atom].x << " " << catom_array[atom].y << " " << catom_array[atom].z << "\tatomj: ";
                  //std::cout << catom_array[natom].x << " " << catom_array[natom].y << " " << catom_array[natom].z << std::endl;

                  // increment 1D counter
                  counter++;
               }
            }
            // push back final start index only if atom is a surface atom
            // [si x x x x x ][ei/si x x x x x][ei/si x x x x x]
            atoms::nearest_neighbour_list_ei[atom]=counter;
         }
         //std::cout << std::endl;
      }

   } // end of surface anisotropy initialisation
   //-------------------------------------------------------------------------------------------------------------

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

   //----------------------------------------------------------------------------
   // Optionally calculate random atomic anisotropy axes for spherical harmonics
   //----------------------------------------------------------------------------
   if(sim::spherical_harmonics && sim::random_anisotropy){
     	// Resize arrays
     	atoms::uniaxial_anisotropy_vector_x.resize(atoms::num_atoms,0.0);
     	atoms::uniaxial_anisotropy_vector_y.resize(atoms::num_atoms,0.0);
     	atoms::uniaxial_anisotropy_vector_z.resize(atoms::num_atoms,0.0);

     	std::vector<double> grain_anisotropy_directions(0);
     	// check for grain level random anisotropy
     	if(grains::random_anisotropy){
        	// resize array storing grain anisotropy vectors
        	grain_anisotropy_directions.resize(3*grains::num_grains);

			// calculate anisotropy directions for all grains on root process
  			if(vmpi::my_rank == 0){
  	  			for(int g=0; g<grains::num_grains; g++){

            	double x = mtrandom::gaussian();
            	double y = mtrandom::gaussian();
            	double z = mtrandom::gaussian();

            	// Calculate vector length
            	const double r = 1.0/sqrt (x*x + y*y + z*z);

            	grain_anisotropy_directions[3*g + 0] = x*r;
            	grain_anisotropy_directions[3*g + 1] = y*r;
            	grain_anisotropy_directions[3*g + 2] = z*r;

				}
			}

			#ifdef MPICF
				// Broadcast calculated anisotropy directions to all nodes
				MPI_Bcast(&grain_anisotropy_directions[0], grain_anisotropy_directions.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#endif

		}

		// Unroll anisotropy directions
		for(int atom = 0; atom < atoms::num_atoms; ++atom){

			// Determine material number
        	int imaterial = atoms::type_array[atom];

		  	// Calculate random anisotropy directions on unit sphere
		  	if(mp::material[imaterial].random_anisotropy){

			  	double x = mtrandom::gaussian();
			  	double y = mtrandom::gaussian();
			  	double z = mtrandom::gaussian();

			  	// Calculate vector length
			  	const double r = 1.0/sqrt (x*x + y*y + z*z);

     			// Save direction
     			atoms::uniaxial_anisotropy_vector_x[atom] = x*r;
     			atoms::uniaxial_anisotropy_vector_y[atom] = y*r;
     			atoms::uniaxial_anisotropy_vector_z[atom] = z*r;

			}
        	// If random grain anisotropy defined set local anisotropy to the grain anisotropy direction
        	else if(mp::material[imaterial].random_grain_anisotropy){

          	const int grain = atoms::grain_array[atom];

          	atoms::uniaxial_anisotropy_vector_x[atom] = grain_anisotropy_directions[3*grain + 0];
          	atoms::uniaxial_anisotropy_vector_y[atom] = grain_anisotropy_directions[3*grain + 1];
          	atoms::uniaxial_anisotropy_vector_z[atom] = grain_anisotropy_directions[3*grain + 2];

			}
			// Otherwise unroll anisotropy directions
			else{
	   		atoms::uniaxial_anisotropy_vector_x[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[0];
	   		atoms::uniaxial_anisotropy_vector_y[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[1];
	   		atoms::uniaxial_anisotropy_vector_z[atom] = mp::material[imaterial].UniaxialAnisotropyUnitVector[2];
			}

		}
	}

   return EXIT_SUCCESS;

}

} // End of cs namespace
