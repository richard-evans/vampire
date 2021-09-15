//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2018. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef HIERARCHICAL_INTERNAL_H_
#define HIERARCHICAL_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the hierarchical module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "hierarchical.hpp"

// hierarchical module headers
#include "internal.hpp"
#include "vector"

namespace hierarchical{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern int num_levels;
      //create arrays for data storage
      extern std::vector < double > cell_positions;
      extern std::vector < double > cell_positions_mom;
      extern std::vector < double > cell_dimensions;
      extern std::vector < int > cells_level_start_index;
      extern std::vector < int > cells_level_end_index;
      extern std::vector < int > cells_in_cells;
      extern std::vector < int > interaction_range;
      extern std::vector < int > cells_in_cells_start_index;
      extern std::vector < int > cells_in_cells_end_index;
      extern std::vector <int> interaction_list;
      extern std::vector < int > interaction_list_start_index;
      extern std::vector < int > interaction_list_end_index;
      extern std::vector < int > num_atoms_in_cell;

      extern std::vector < int > proc_cell_index_array1D;
      extern std::vector < int > cells_num_atoms_in_cell;

      extern int total_num_cells;
      extern int num_zero_level_cells;

      extern std::vector<double> mag_array_x; /// arrays to store cells magnetisation
      extern std::vector<double> mag_array_y;
      extern std::vector<double> mag_array_z;

      extern std::vector<double> rij_tensor_xx;
      extern std::vector<double> rij_tensor_xy;
      extern std::vector<double> rij_tensor_xz;

      extern std::vector<double> rij_tensor_yy;
      extern std::vector<double> rij_tensor_yz;
      extern std::vector<double> rij_tensor_zz;


      extern double av_cell_size;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      std::vector < std::vector < double> > calculate_corners(double x, double y, double z, double cell_size_x, double cell_size_y, double cell_size_z);

      // function to update hierarchical cell magnetizations
      void calculate_hierarchical_magnetisation(std::vector <double>& x_spin_array, // atomic spin directions
                                                std::vector <double>& y_spin_array,
                                                std::vector <double>& z_spin_array,
                                                std::vector <double>& m_spin_array, // atomic spin moment
                                                std::vector < bool >& magnetic);    // is magnetic);

      // new version of inter tensor method
      void calc_inter(const int celli,                                                // global ID of cell i
                      const int cellj,                                                // global ID of cell i
                      const int interaction_no,                                       // interaction number
                      const double cutoff,                                            // cutoff range for dipole tensor construction (Angstroms)
                      const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                      const std::vector<double>& cells_pos_and_mom_array,             // array of positions and cell moments
                      const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                      const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                     );

      void calc_intra(const int celli,                                                // global ID of cell i
                      const int cellj,                                                // global ID of cell i
                      const int interaction_no,                                       // interaction number
                      const std::vector<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                      const std::vector<int>& list_of_cells_with_atoms,               // list of cells to access atoms
                      const std::vector< std::vector<double> >& atoms_in_cells_array  // output array of positions and moments of atoms in cells
                     );

      void calc_tensor(int cells_num_local_cells,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_x,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_y,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_z);

   } // end of internal namespace

} // end of hierarchical namespace

#endif //HIERARCHICAL_INTERNAL_H_
