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
#ifndef ATOMS_H_
#define ATOMS_H_

#include <string>
#include <vector>

#include "exchange.hpp"

// unit vector type
class uvec_t{
   public:
      double x; // x hat
      double y; // y hat
      double z; // z hat
      double r; // length of vector

      // unit dot product
      double dot(uvec_t a, uvec_t b){
         return (a.x * b.x + a.y * b.y + a.z * b.z);
      }

      // dot product including length
      double rdot(uvec_t a, uvec_t b){
         return (a.r * b.r * (a.x * b.x + a.y * b.y + a.z * b.z));
      }

      // unit cross product
      uvec_t cross(uvec_t a, uvec_t b){
         uvec_t tmp;
         tmp.x = a.y * b.z - a.z * b.y;
         tmp.y = a.z * b.x - a.x * b.z;
         tmp.z = a.x * b.y - a.y * b.x;
         tmp.r = a.r * b.r;

         return tmp;
      }

};

//======================================================================
//                       Global Atomistic Variables
//======================================================================
namespace atoms
{
	//--------------------------
	// Single Variables
	//--------------------------
	extern int num_atoms;			/// Number of atoms in simulation
	extern int num_neighbours;	   	/// Maximum number of neighbours for Hamiltonian/Lattice
	extern int total_num_neighbours;/// Total number of neighbours for system
   extern uint64_t num_non_magnetic_atoms; // Number of non-magnetic atoms not to be simulated

	//--------------------------
	// Array Variables
	//--------------------------

	extern std::vector <double> x_coord_array;
	extern std::vector <double> y_coord_array;
	extern std::vector <double> z_coord_array;
	extern std::vector <int> neighbour_list_array;
	extern std::vector <int> neighbour_interaction_type_array;
	extern std::vector <int> neighbour_list_start_index;
	extern std::vector <int> neighbour_list_end_index;
	extern std::vector <int> type_array;
	extern std::vector <int> category_array;
	extern std::vector <int> grain_array;
	extern std::vector <int> cell_array;

	extern std::vector <double> x_spin_array;
	extern std::vector <double> y_spin_array;
	extern std::vector <double> z_spin_array;
   extern std::vector <double> m_spin_array; /// Array of atomic spin moments

	extern std::vector <double> x_total_spin_field_array;		/// Total spin dependent fields
	extern std::vector <double> y_total_spin_field_array;		/// Total spin dependent fields
	extern std::vector <double> z_total_spin_field_array;		/// Total spin dependent fields
	extern std::vector <double> x_total_external_field_array;	/// Total external fields
	extern std::vector <double> y_total_external_field_array;	/// Total external fields
	extern std::vector <double> z_total_external_field_array;	/// Total external fields

	extern std::vector <zval_t> i_exchange_list;
	extern std::vector <zvec_t> v_exchange_list;
	extern std::vector <zten_t> t_exchange_list;

   extern std::vector <bool> surface_array; // flag to identify atom as surface
   extern std::vector <bool> magnetic; // flag to identify atom as being magnetic

   extern std::vector <uvec_t> neighbour_eij_array; // unrolled list of eij unit vectors between neighbouring atoms

}


#endif /*ATOMS_H_*/
