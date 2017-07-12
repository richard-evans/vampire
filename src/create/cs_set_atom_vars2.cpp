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

#include "anisotropy.hpp"
#include "atoms.hpp"
#include "cells.hpp"
#include "create.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "exchange.hpp"
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

   //---------------------------------------------------------------------------
   // Identify surface atoms and initialise anisotropy data
   //---------------------------------------------------------------------------
   anisotropy::identify_surface_atoms(catom_array, cneighbourlist);

	//===========================================================
	// Create 1-D neighbourlist
	//===========================================================

	zlog << zTs() << "Memory required for creation of 1D neighbour list on rank " << vmpi::my_rank << ": ";
	zlog << (2.0*double(atoms::num_atoms)+2.0*double(atoms::total_num_neighbours))*8.0/1.0e6 << " MB RAM"<< std::endl;

   //-------------------------------------------------
	//	Initialise exchange calculation
	//-------------------------------------------------
   exchange::initialize(cneighbourlist);

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
