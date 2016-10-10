//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef DIPOLE_H_
#define DIPOLE_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "dipole.hpp"
#include "vmpi.hpp"
#include "material.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for dipole module
//--------------------------------------------------------------------------------
namespace dipole{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   extern int update_rate; /// timesteps between updates
   extern bool activated;

   extern std::vector<double> cells_field_array_x; /// arrays to store cells field
   extern std::vector<double> cells_field_array_y;
   extern std::vector<double> cells_field_array_z;
   extern std::vector < double > atom_dipolar_field_array_x;
   extern std::vector < double > atom_dipolar_field_array_y;
   extern std::vector < double > atom_dipolar_field_array_z;

   extern double cutoff;

   //-----------------------------------------------------------------------------
   // Function to unroll cells dipolar field into atomic field
   //-----------------------------------------------------------------------------
   void calculate_field();

   ////-----------------------------------------------------------------------------
   //// Function to send receive cells and atoms data to other cpus
   ////-----------------------------------------------------------------------------
   //void send_receive_data(int num_local_atoms, std::vector<double>& atom_coords_x, std::vector<double>& atom_coords_y, std::vector<double>& atom_coords_z);
   ////-----------------------------------------------------------------------------
   //// Function to locate cells cpus
   ////-----------------------------------------------------------------------------
   //void locate_cells_on_cpu(int cpu,int local_cell,double x,double y,double z,std::vector<int>& proc_index_array1D,std::vector<int>& proc_cell_index_array1D,std::vector< std::vector<int> >& proc_cell_index_array,std::vector<double>& minimax);
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
                           int cells_num_cells
                           );
   //-----------------------------------------------------------------------------
   // Function to send receive atoms data to other cpus
   //-----------------------------------------------------------------------------
   //int send_recv_atoms_data(std::vector<int>& proc_cell_index_array1D,std::vector<int>& cell_id_array,std::vector<double>& atom_pos_x,std::vector<double>& atom_pos_y,std::vector<double>& atom_pos_z,std::vector<int>& atom_type_array,int cell,int cpu_send,int cpu_recv,MPI::Status status,MPI::Request send_request,MPI::Request recv_request,int num_local_atoms);
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
                           double cells_macro_cell_size
                           );
   /*------------------------------------------------*/
   /*Function to sort cells/atoms data after sharing */
   /*------------------------------------------------*/
   int sort_data(std::vector<int>& proc_cell_index_array1D,
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

   //-----------------------------------------------------------------------------
   // Function to initialise dipole module
   //-----------------------------------------------------------------------------
   void initialize(const int cells_num_atoms_in_unit_cell,
                   int cells_num_cells, /// number of macrocells
                   int cells_num_local_cells, /// number of local macrocells
                   const double cells_macro_cell_size,
                   std::vector <int>& cells_local_cell_array,
                   std::vector <int>& cells_num_atoms_in_cell, /// number of atoms in each cell
                   std::vector < std::vector <int> >& cells_index_atoms_array,
                   const std::vector<double>& cells_volume_array,
                   /*const std::vector<double>& cells_cell_coords_array_x, /// arrays to store cells positions
                   const std::vector<double>& cells_cell_coords_array_y,
                   const std::vector<double>& cells_cell_coords_array_z,*/

                   std::vector<double>& cells_pos_and_mom_array,  // array to store positions and moment of cells

                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                   const std::vector<int>& atom_type_array,
                   const std::vector<int>& atom_cell_id_array,

                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,

                   const int num_atoms
   );

   //---------------------------------------------------------------------------
   // Function to process input file parameters for dipole module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, std::vector<mp::materials_t>& read_material);

} // end of dipole namespace

#endif //DIPOLE_H_
