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
#include <vector>

namespace stats
//==========================================================
// Namespace statistics
//==========================================================
{
	extern int num_atoms;				// Number of atoms for statistic purposes
	extern double inv_num_atoms;	//1.0/num_atoms
	extern double max_moment;		// Total Maximum moment

	extern double total_mag_actual[3];	///< Actual magnetisation components
	extern double total_mag_m_actual;	///< Actual magnitude of total magnetisation
	extern double total_mean_mag_m_actual;	///< Actual magnitude of total magnetisation
	
	extern double total_mag_norm[3];	///< Normalised magnetisation components
	extern double total_mag_m_norm;	///< Normalised magnitude of total magnetisation
	extern double total_mean_mag_m_norm;	///< Normalised magnitude of total magnetisation

	extern double data_counter;		// number of data points for averaging

	// Member Functions
	extern int mag_m();
	extern void mag_m_reset();
	extern double max_torque();
	
	extern std::vector <double> sublattice_mx_array;
	extern std::vector <double> sublattice_my_array;
	extern std::vector <double> sublattice_mz_array;
	extern std::vector <double> sublattice_magm_array;
	extern std::vector <double> sublattice_mean_magm_array;
	extern std::vector <double> sublattice_mom_array;
	extern std::vector <int> sublattice_nm_array;
	
	extern bool calculate_torque;
	extern double total_system_torque[3];
	extern double total_mean_system_torque[3];
	
	extern std::vector <double> sublattice_mean_torque_x_array;
	extern std::vector <double> sublattice_mean_torque_y_array;
	extern std::vector <double> sublattice_mean_torque_z_array;

	extern double torque_data_counter;

   extern double mean_susceptibility[3];
   extern double mean_susceptibility_squared[3];
   extern bool calculate_susceptibility;

   extern bool calculate_energy;

   // Statistics energy types
   enum energy_t { all=0, exchange=1, anisotropy=2, cubic_anisotropy=3, surface_anisotropy=4,applied_field=5, magnetostatic=6, second_order_anisotropy=7 };

   // Statistics types
   enum stat_t { total=0, mean=1};

   // Statistics output functions
   extern void output_energy(std::ostream&, enum energy_t, enum stat_t);

}
