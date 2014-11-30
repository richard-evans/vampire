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
#ifndef STATS_H_
#define STATS_H_
#include <vector>
#include <string>

namespace stats
//==========================================================
// Namespace statistics
//==========================================================
{
	extern int num_atoms;				/// Number of atoms for statistic purposes
	extern double inv_num_atoms;	///1.0/num_atoms
	extern double max_moment;		/// Total Maximum moment
	extern double data_counter;		/// number of data points for averaging

	// Member Functions
	extern int mag_m();
	extern void mag_m_reset();
	extern double max_torque();

	extern bool calculate_torque;
	extern double total_system_torque[3];
	extern double total_mean_system_torque[3];

	extern std::vector <double> sublattice_mean_torque_x_array;
	extern std::vector <double> sublattice_mean_torque_y_array;
	extern std::vector <double> sublattice_mean_torque_z_array;

	extern double torque_data_counter;

   extern double mean_susceptibility[4];
   extern double mean_susceptibility_squared[4];
   extern bool calculate_susceptibility;

   extern bool calculate_energy;

   /// Statistics energy types
   enum energy_t { all=0, exchange=1, anisotropy=2, cubic_anisotropy=3, surface_anisotropy=4,applied_field=5, magnetostatic=6, second_order_anisotropy=7 };

   /// Statistics types
   enum stat_t { total=0, mean=1};

   /// Statistics output functions
   extern void output_energy(std::ostream&, enum energy_t, enum stat_t);

   //-------------------------------------------------
   // New statistics module functions and variables
   //-------------------------------------------------

   // Control functions
   void initialize(const int num_atoms, const int num_materials, const std::vector<double>& magnetic_moment_array, 
                   const std::vector<int>& material_type_array, const std::vector<int>& height_category_array);
   void update(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& mm);
   void reset();

   // Statistics control flags (to be moved internally when long-awaited refactoring of vio is done)
   extern bool calculate_system_magnetization;
   extern bool calculate_material_magnetization;
   extern bool calculate_height_magnetization;
   extern bool calculate_material_height_magnetization;

   //----------------------------------
   // Magnetization Class definition
   //----------------------------------
   class magnetization_statistic_t{

      public:
         //magnetization_statistic_t (const int in_mask_size, std::vector<int> in_mask);
         magnetization_statistic_t ();
         void set_mask(const int mask_size, std::vector<int> inmask, const std::vector<double>& mm);
         void calculate_magnetization(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& mm);
         void reset_magnetization_averages();
         const std::vector<double>& get_magnetization();
         std::string output_magnetization();
         std::string output_normalized_magnetization();
         std::string output_normalized_magnetization_length();
         std::string output_normalized_mean_magnetization();
         std::string output_normalized_mean_magnetization_length();
         std::string output_normalized_magnetization_dot_product(const std::vector<double>& vec);

      private:
         bool is_initialized;
         int num_atoms;
         int mask_size;
         double mean_counter;
         std::vector<int> mask;
         std::vector<double> magnetization;
         std::vector<double> mean_magnetization;
         std::vector<int> zero_list;
         std::vector<double> saturation;

   };

   // Statistics classes
   extern magnetization_statistic_t system_magnetization;
   extern magnetization_statistic_t material_magnetization;
   extern magnetization_statistic_t height_magnetization;
   extern magnetization_statistic_t material_height_magnetization;


}

#endif /*STATS_H_*/
