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
         atomistic = 4 // new macrocell with tensor including local corrections and nearfield multipole
         //exact        = 4, // atomistic dipole dipole (too slow for anything over 1000 atoms)
      };

      extern solver_t solver;

      extern int update_time; /// last update time

      extern const double prefactor; // 1e-7/1e30

      extern std::vector <std::vector < double > > rij_tensor_xx;
      extern std::vector <std::vector < double > > rij_tensor_xy;
      extern std::vector <std::vector < double > > rij_tensor_xz;

      extern std::vector <std::vector < double > > rij_tensor_yy;
      extern std::vector <std::vector < double > > rij_tensor_yz;
      extern std::vector <std::vector < double > > rij_tensor_zz;

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

      //------------------------------------------------------------------------
      // data structures for atomistic solver
      // (copy of all atom positions and spins on all processors)
      //------------------------------------------------------------------------

      extern int num_local_atoms; // number of local atoms (my processor)
      extern int total_num_atoms; // number of total atoms (all processors)

      // arrays to store atomic coordinates
      extern std::vector <double> cx;
      extern std::vector <double> cy;
      extern std::vector <double> cz;

      // arrays to store atomic spins
      extern std::vector <double> sx;
      extern std::vector <double> sy;
      extern std::vector <double> sz;
      extern std::vector <double> sm;

      // arrays for calculating displacements for parallelisation
      extern std::vector <int> receive_counts;
      extern std::vector <int> receive_displacements;

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
                                    std::vector<double>& cells_volume_array,
                                    std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                    std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                    std::vector<int>& atom_type_array,
                                    std::vector<int>& atom_cell_id_array,
                                    std::vector<double>& atom_coords_x, //atomic coordinates
                                    std::vector<double>& atom_coords_y,
                                    std::vector<double>& atom_coords_z,
                                    int num_atoms);

      void compute_inter_tensor(const double cells_macro_cell_size,
                                const int i,
                                const int j,
                                const int lc,
                                std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                //std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z);

      void compute_intra_tensor(const int i,
                                const int j,
                                const int lc,
                                std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z);

      void initialize_macrocell_solver(const int cells_num_atoms_in_unit_cell,
                                       int cells_num_cells, /// number of macrocells
                                       int cells_num_local_cells, /// number of local macrocells
                                       const double cells_macro_cell_size,
                                       std::vector <int>& cells_local_cell_array,
                                       std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                                       std::vector <int>& cells_num_atoms_in_cell_global, /// number of atoms in each cell
                                       std::vector < std::vector <int> >& cells_index_atoms_array,
                                       std::vector<double>& cells_volume_array,
                                       std::vector<double>& cells_pos_and_mom_array, // array to store positions and moment of cells
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                                       std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                                       std::vector<int>& atom_type_array,
                                       std::vector<int>& atom_cell_id_array,
                                       std::vector<double>& atom_coords_x, //atomic coordinates
                                       std::vector<double>& atom_coords_y,
                                       std::vector<double>& atom_coords_z,
                                       int num_atoms);

      void initialize_atomistic_solver(int num_atoms,                      // number of atoms (only correct in serial)
                                       std::vector<double>& x_coord_array, // atomic corrdinates (angstroms)
                                       std::vector<double>& y_coord_array,
                                       std::vector<double>& z_coord_array,
                                       std::vector<double>& moments_array); // atomistic magnetic moments (bohr magnetons)

      void calculate_atomistic_dipole_field(std::vector<double>& x_spin_array, // atomic spin directions
                                            std::vector<double>& y_spin_array,
                                            std::vector<double>& z_spin_array);

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

      /*--------------------------------------------------------*/
      /*Function to send cells field to be output in cfg file   */
      /*--------------------------------------------------------*/
      int send_cells_demag_factor(std::vector<int>& cells_cell_id_array,
                                 std::vector<double>& N_tensor_array,
                                 int cells_num_local_cells);

   } // end of internal namespace

} // end of dipole namespace

#endif //DIPOLE_INTERNAL_H_
