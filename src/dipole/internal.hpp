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
#include "vmpi.hpp"

// dipole module headers
#include "internal.hpp"

namespace dipole{

   extern int update_rate; /// timesteps between updates

   extern void calculate_field();
   //extern void send_receive_data(int num_local_atoms,std::vector<double>& atom_coords_x, std::vector<double>& atom_coords_y, std::vector<double>& atom_coords_z);
   //extern void locate_cells_on_cpu(int cpu,int local_cell,double x,double y,double z,std::vector<int>& proc_index_array1D,std::vector<int>& proc_cell_index_array1D,std::vector< std::vector<int> >& proc_cell_index_array,std::vector<double>& minimax);
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

   //extern int send_recv_atoms_data(std::vector<int>& proc_cell_index_array2D,std::vector<int>& cell_id_array,std::vector<double>& atom_pos_x,std::vector<double>& atom_pos_y,std::vector<double>& atom_pos_z,std::vector<int>& atom_type_array,int cell,int cpu_send,int cpu_recv,MPI::Status status,MPI::Request send_request,MPI::Request recv_request,int num_local_atoms);
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
                                 int cells_macro_cell_size
                                 );

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool initialised;
      extern bool enabled;

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

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      //void write_macrocell_data();
      extern void update_field();

   } // end of internal namespace

} // end of dipole namespace

#endif //DIPOLE_INTERNAL_H_
