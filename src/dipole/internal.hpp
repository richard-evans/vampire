//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef DIPOLE_INTERNAL_H_
#define DIPOLE_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the dipole module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "dipole.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool initialised;

      // enumerated list of different dipole solvers
      enum solver_t{
         macrocell    = 0, // original bare macrocell method (cheap but inaccurate)
         tensor       = 1, // new macrocell with tensor including local corrections
         //multipole    = 2, // bare macrocell but with multipole expansion
         //hierarchical = 3, // new macrocell with tensor including local corrections and nearfield multipole
         //exact        = 4, // atomistic dipole dipole (too slow for anything over 1000 atoms)
      };

      extern solver_t solver;

      extern int update_time; /// last update time

      extern const double prefactor; // 1e-7/1e30

      extern std::vector <std::vector < double > > rij_inter_xx;
      extern std::vector <std::vector < double > > rij_inter_xy;
      extern std::vector <std::vector < double > > rij_inter_xz;

      extern std::vector <std::vector < double > > rij_inter_yy;
      extern std::vector <std::vector < double > > rij_inter_yz;
      extern std::vector <std::vector < double > > rij_inter_zz;

      extern std::vector <std::vector < double > > rij_intra_xx;
      extern std::vector <std::vector < double > > rij_intra_xy;
      extern std::vector <std::vector < double > > rij_intra_xz;

      extern std::vector <std::vector < double > > rij_intra_yy;
      extern std::vector <std::vector < double > > rij_intra_yz;
      extern std::vector <std::vector < double > > rij_intra_zz;

      extern int num_atoms;
      extern std::vector < int > atom_type_array;
      extern std::vector < int > atom_cell_id_array;

      extern int cells_num_cells;
      extern int cells_num_local_cells;
      extern std::vector <int>  cells_local_cell_array;
      extern std::vector <int>  cells_num_atoms_in_cell;
      extern std::vector < double > cells_volume_array;

      extern std::vector<double> cells_pos_and_mom_array;
      extern std::vector < int > proc_cell_index_array1D;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      //void write_macrocell_data();
      extern void update_field();

      void allocate_memory(const int cells_num_local_cells, const int cells_num_cells);

      void initialize_tensor_solver(const int cells_num_atoms_in_unit_cell,
                                    int cells_num_cells, /// number of macrocells
                                    int cells_num_local_cells, /// number of local macrocells
                                    const double cells_macro_cell_size,
                                    std::vector <int>& cells_local_cell_array,
                                    std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                    std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                    std::vector < std::vector <int> >& cells_index_atoms_array,
                                    const std::vector<double>& cells_volume_array,
                                    std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                    const std::vector<int>& atom_type_array,
                                    const std::vector<int>& atom_cell_id_array,
                                    const std::vector<double>& atom_coords_x, //atomic coordinates
                                    const std::vector<double>& atom_coords_y,
                                    const std::vector<double>& atom_coords_z,
                                    const int num_atoms);

      void initialize_macrocell_solver();

      //-----------------------------------------------------------------------------
      // Function to send receive cells data to other cpus
      //-----------------------------------------------------------------------------
      int send_recv_cells_data(std::vector<int>& proc_cell_index_array1D,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                               std::vector< std::vector <int> >& cells_index_atoms_array,
                               std::vector<double>& cells_pos_and_mom_array,
                               std::vector<int>& cells_num_atoms_in_cell,
                               std::vector<int>& cells_cell_id_array,
                               std::vector<int>& cells_local_cell_array,
                               int cells_num_local_cells,
                               int cells_num_cells);

      //-----------------------------------------------------------------------------
      // Function to send receive atoms data to other cpus
      //-----------------------------------------------------------------------------
      int send_recv_atoms_data(std::vector<int>& proc_cell_index_array2D,
                               std::vector<int>& cell_id_array,
                               std::vector<int>& cells_local_cell_array,
                               std::vector<double>& atom_pos_x,
                               std::vector<double>& atom_pos_y,
                               std::vector<double>& atom_pos_z,
                               std::vector<int>& atom_type_array,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                               std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                               std::vector< std::vector <int> >& cells_index_atoms_array,
                               std::vector<double>& cells_pos_and_mom_array,
                               std::vector<int>& cells_num_atoms_in_cell,
                               int cells_num_local_cells,
                               int cells_num_cells,
                               double cells_macro_cell_size);

      //----------------------------------------------------------------
      //Function to sort cells/atoms data after sharing
      //----------------------------------------------------------------
      int sort_data(std::vector<int>& proc_cell_index_array1D,
                  std::vector<int>& cells_cell_id_array,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                  std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                  std::vector< std::vector <int> >& cells_index_atoms_array,
                  std::vector<double>& cells_pos_and_mom_array,
                  std::vector<int>& cells_num_atoms_in_cell,
                  int cells_num_local_cells,
                  int cells_num_cells);

   } // end of internal namespace

} // end of dipole namespace

#endif //DIPOLE_INTERNAL_H_
