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

   extern int update_rate; /// timesteps between updates
   extern void calculate_field();
   extern int send_recv_cells_data(std::vector<int>& proc_cell_index_array1D,
                                 std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                 std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                 std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                 std::vector< std::vector <int> >& cells_index_atoms_array,
                                 std::vector<double>& cells_pos_and_mom_array,
                                 std::vector<int>& cells_num_atoms_in_cell,
                                 std::vector<int>& cells_cell_id_array,
                                 std::vector<int>& cells_local_cell_array,
                                 int cells_num_local_cells,
                                 int cells_num_cells
                                 );

   extern int send_recv_atoms_data(std::vector<int>& proc_cell_index_array2D,
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
                                 double cells_macro_cell_size
                                 );

   extern int sort_data(std::vector<int>& proc_cell_index_array1D,
                        std::vector<int>& cells_cell_id_array,
                        std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_x,
                        std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_y,
                        std::vector< std::vector <double> >& cells_atom_in_cell_coords_array_z,
                        std::vector< std::vector <int> >& cells_index_atoms_array,
                        std::vector<double>& cells_pos_and_mom_array,
                        std::vector<int>& cells_num_atoms_in_cell,
                        int cells_num_local_cells,
                        int cells_num_cells
				  		   );

   extern int send_cells_field(std::vector<int>& cells_cell_id_array,
                              std::vector<double>& dipole_cells_field_array_x,
                              std::vector<double>& dipole_cells_field_array_y,
                              std::vector<double>& dipole_cells_field_array_z,
                              std::vector<double>& dipole_cells_mu0Hd_field_array_x, // mu_0*Hd-field
                              std::vector<double>& dipole_cells_mu0Hd_field_array_y,
                              std::vector<double>& dipole_cells_mu0Hd_field_array_z,
                              std::vector<double>& cells_volume_array,
                              int cells_num_local_cells
                  );

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool initialised;

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

   } // end of internal namespace

} // end of dipole namespace

#endif //DIPOLE_INTERNAL_H_
