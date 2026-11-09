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
#include <string>

// Vampire headers
#include "hierarchical.hpp"
#include "vmpi.hpp"

// hierarchical module headers
#include "vector"

namespace hierarchical{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern int num_levels;
      // arrays for data storage; the shared arrays are held once per shared memory region
      // (cell geometry is needed only during initialisation and freed at its end)
      extern vmpi::shared_array_t < double > cell_positions;
      extern vmpi::shared_array_t < double > cell_positions_mom;
      extern vmpi::shared_array_t < double > cell_dimensions;
      extern std::vector < int > cells_level_start_index;
      extern std::vector < int > cells_level_end_index;
      extern vmpi::shared_array_t < int > cells_in_cells;
      extern std::vector < int > interaction_range;
      extern vmpi::shared_array_t < int > cells_in_cells_start_index;
      extern vmpi::shared_array_t < int > cells_in_cells_end_index;
      extern std::vector <int> interaction_list;
      extern std::vector < int > interaction_list_start_index;
      extern std::vector < int > interaction_list_end_index;
      extern std::vector < int > upward_pass_cells; // coarse cells (level >= 1, by level) holding this process's cells

      extern std::vector < int > proc_cell_index_array1D;
      extern std::vector < int > cells_num_atoms_in_cell;

      extern int total_num_cells;
      extern int num_zero_level_cells;

      extern vmpi::shared_array_t<double> mag_array_x; /// cell magnetisations (all processes), one copy per region
      extern vmpi::shared_array_t<double> mag_array_y;
      extern vmpi::shared_array_t<double> mag_array_z;
      extern std::vector<int> contributing_cells;      // cells (all levels, in order) with a contribution from this process
      extern std::vector<double> partial_mag_array_x;  /// this process's contribution to the moment of each of contributing_cells
      extern std::vector<double> partial_mag_array_y;  /// (compact: indexed like contributing_cells, not by cell id)
      extern std::vector<double> partial_mag_array_z;
      extern std::vector<int> atom_contributing_index;          // index in contributing_cells of each local atom's cell (-1 if not magnetic)
      extern std::vector<int> micromagnetic_contributing_index; // index in contributing_cells of each micromagnetic cell of this process
      extern std::vector<int> upward_pass_index;       // index in contributing_cells of each of upward_pass_cells
      extern std::vector<int> upward_children_start;   // children of upward_pass_cells[i] in contributing_cells:
      extern std::vector<int> upward_children;         //    upward_children[upward_children_start[i] .. upward_children_start[i+1])
      extern std::vector<double> cell_field_x;   // field of each local cell (indexed like cells::cell_id_array)
      extern std::vector<double> cell_field_y;
      extern std::vector<double> cell_field_z;
      extern std::vector<int> atom_local_cell;   // local cell index of each local magnetic atom's cell (-1 otherwise)
      extern bool full_cell_fields; // also fill dipole::cells_field_array_* (sized for all cells)
      extern std::vector<double> mag_reduction_values; // contributions of this process, x then y then z
      extern vmpi::shared_sum_t mag_reduction;         // sum of the contributions over processes

      extern std::vector<double> rij_tensor_xx;
      extern std::vector<double> rij_tensor_xy;
      extern std::vector<double> rij_tensor_xz;

      extern std::vector<double> rij_tensor_yy;
      extern std::vector<double> rij_tensor_yz;
      extern std::vector<double> rij_tensor_zz;


      extern double av_cell_size;

      extern bool output_p2p_diagnostic; // log the cost of a point-to-point exchange of cell moments at initialisation

      // how the cell moments are summed between shared memory regions (hierarchical:moment-exchange)
      enum class moment_exchange_t { all_reduce, point_to_point };
      extern moment_exchange_t moment_exchange;
      extern int moment_exchange_level; // highest level exchanged point to point (hierarchical:moment-exchange-level)

      // phases of the field update, timed separately for the scaling study (papers/dipole)
      enum phase_t { upward_pass = 0, moment_reduction = 1, field_loop = 2, field_distribution = 3, num_phases = 4 };
      extern std::vector<double> phase_time; // time (s) of each phase of the last update on this process

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      std::vector < std::vector < double> > calculate_corners(double x, double y, double z, double cell_size_x, double cell_size_y, double cell_size_z);

      // function to log the cost of a point-to-point exchange of cell moments
      void write_p2p_diagnostic();

      // function to set up the point-to-point exchange of cell moments between regions (exchange.cpp)
      void initialise_moment_exchange();
      void exchange_moments();                   // sum of the cell moments, point to point (in place of mag_reduction.sum)
      std::string moment_exchange_phase_times(); // phase times of the last exchange for the DIPOLE_SCALING line

      // function to update hierarchical cell magnetizations
      void calculate_hierarchical_magnetisation(std::vector <double>& x_spin_array, // atomic spin directions
                                                std::vector <double>& y_spin_array,
                                                std::vector <double>& z_spin_array,
                                                std::vector <double>& m_spin_array, // atomic spin moment
                                                std::vector < bool >& magnetic);    // is magnetic);

      // Table of the atoms of each level 0 cell, one per shared memory region (region_atoms.cpp): the magnetic
      // atoms of the region's processes and of other regions' cells within the atomistic cutoff, sorted by cell,
      // source process and index. The atoms of the k-th cell held, cell_ids[k], are a = cell_starts[k] ..
      // cell_starts[k+1]-1, with coordinates coords[3a+0..2].
      struct region_atoms_t{
         vmpi::shared_array_t<double> coords;
         vmpi::shared_array_t<int> cell_ids;    // ascending
         vmpi::shared_array_t<int> cell_starts; // one more than cell_ids
         void build(const double cutoff, const std::vector<double>& atom_coords_x, const std::vector<double>& atom_coords_y,
                    const std::vector<double>& atom_coords_z); // collective
         int find(const int cell) const;        // index k of cell, or -1
         void free();                           // collective
      };
      extern region_atoms_t region_atoms;

      // new version of inter tensor method
      void calc_inter(const int celli,                                                // global ID of cell i
                      const int cellj,                                                // global ID of cell i
                      const int interaction_no,                                       // interaction number
                      const double cutoff,                                            // cutoff range for dipole tensor construction (Angstroms)
                      const vmpi::shared_array_t<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                      const region_atoms_t& atoms_table                               // atoms of the cells near this process's cells
                     );

      void calc_intra(const int celli,                                                // global ID of cell i
                      const int cellj,                                                // global ID of cell i
                      const int interaction_no,                                       // interaction number
                      const vmpi::shared_array_t<int>& global_atoms_in_cell_count,             // number of atoms in each cell (all CPUs)
                      const region_atoms_t& atoms_table                               // atoms of the cells near this process's cells
                     );

      void calc_tensor(int cells_num_local_cells,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_x,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_y,
                       std::vector < std::vector < double > >& cells_atom_in_cell_coords_array_z);

   } // end of internal namespace

} // end of hierarchical namespace

#endif //HIERARCHICAL_INTERNAL_H_
