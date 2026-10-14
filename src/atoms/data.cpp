//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#include "atoms.hpp"

#include <vector>
#include <cstdint>

//==========================================================
// Namespace atom variables
//==========================================================
namespace atoms{
   //--------------------------
   // Single Variables
   //--------------------------
   int num_atoms = 0;         /// Number of atoms in simulation
   int num_neighbours = 0;       /// Maximum number of neighbours for Hamiltonian/Lattice
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
   std::vector <double> mass_spin_array(0);

   std::vector <double> x_total_spin_field_array(0);     /// Total spin dependent fields
   std::vector <double> y_total_spin_field_array(0);     /// Total spin dependent fields
   std::vector <double> z_total_spin_field_array(0);     /// Total spin dependent fields

   std::vector <double> x_total_spin_forces_array(0);    /// Total spin dependent forces
   std::vector <double> y_total_spin_forces_array(0);    /// Total spin dependent forces
   std::vector <double> z_total_spin_forces_array(0);    /// Total spin dependent forces

   std::vector <double> x_velo_array(0);     ///  velocities
   std::vector <double> y_velo_array(0);     ///
   std::vector <double> z_velo_array(0);     ///

   std::vector <double> x_total_external_field_array(0); /// Total external fields
   std::vector <double> y_total_external_field_array(0); /// Total external fields
   std::vector <double> z_total_external_field_array(0); /// Total external fields

   std::vector <double> x_thermal_field_array(0);  /// Total external fields
   std::vector <double> y_thermal_field_array(0);  /// Total external fields
   std::vector <double> z_thermal_field_array(0);  /// Total external fields


   std::vector <zval_t> i_exchange_list(0);
   std::vector <zvec_t> v_exchange_list(0);
   std::vector <zten_t> t_exchange_list(0);

   std::vector <bool> surface_array(0); // flag to identify atom as surface
   std::vector <bool> magnetic(0); // flag to identify atom as being magnetic

   std::vector <uvec_t> neighbour_eij_array; // unrolled list of eij unit vectors between neighbouring atoms

}
