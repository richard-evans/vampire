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
#include "atoms.hpp"

#include <vector>

//==========================================================
// Namespace atom variables
//==========================================================
namespace atoms{
	//--------------------------
	// Single Variables
	//--------------------------
   int num_atoms = 0;			/// Number of atoms in simulation
   int num_neighbours = 0;	   	/// Maximum number of neighbours for Hamiltonian/Lattice
   int total_num_neighbours = 0;
   uint64_t num_non_magnetic_atoms = 0; // Number of non-magnetic atoms not to be simulated

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
    std::vector <double> m_spin_array(0);

	std::vector <double> x_total_spin_field_array(0);		/// Total spin dependent fields
	std::vector <double> y_total_spin_field_array(0);		/// Total spin dependent fields
	std::vector <double> z_total_spin_field_array(0);		/// Total spin dependent fields
	std::vector <double> x_total_external_field_array(0);	/// Total external fields
	std::vector <double> y_total_external_field_array(0);	/// Total external fields
	std::vector <double> z_total_external_field_array(0);	/// Total external fields

	std::vector <zval_t> i_exchange_list(0);
	std::vector <zvec_t> v_exchange_list(0);
	std::vector <zten_t> t_exchange_list(0);

   std::vector <bool> surface_array(0); // flag to identify atom as surface
   std::vector <bool> magnetic(0); // flag to identify atom as being magnetic

   std::vector <uvec_t> neighbour_eij_array; // unrolled list of eij unit vectors between neighbouring atoms

}
