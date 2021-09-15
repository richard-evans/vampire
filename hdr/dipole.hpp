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

//--------------------------------------------------------------------------------
// Namespace for variables and functions for dipole module
//--------------------------------------------------------------------------------
namespace dipole{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   extern bool activated;

   extern int update_rate; /// timesteps between updates
   extern int update_time; /// last update time

   extern std::vector<double> cells_field_array_x; /// arrays to store cells field
   extern std::vector<double> cells_field_array_y;
   extern std::vector<double> cells_field_array_z;
   extern std::vector<double> atom_dipolar_field_array_x;
   extern std::vector<double> atom_dipolar_field_array_y;
   extern std::vector<double> atom_dipolar_field_array_z;

   extern std::vector<int> atomistic_dd_neighbourlist;
   extern std::vector<int> atomistic_dd_neighbourlist_start;
   extern std::vector<int> atomistic_dd_neighbourlist_end;

   extern std::vector<double> cells_mu0Hd_field_array_x;             /// arrays to store cells mu_0*Hdemag-field
   extern std::vector<double> cells_mu0Hd_field_array_y;
   extern std::vector<double> cells_mu0Hd_field_array_z;
   extern std::vector<double> atom_mu0demag_field_array_x;        /// arrays to store atoms mu_0*Hdemag-field
   extern std::vector<double> atom_mu0demag_field_array_y;
   extern std::vector<double> atom_mu0demag_field_array_z;

   extern std::vector <int> dipole_cells_num_atoms_in_cell;             /// Array to store number of atoms in cells that will be needed to print the cell config file

   extern double cutoff;
   extern double atomistic_cutoff;
   extern bool atomsitic_tensor_enabled;

   //-----------------------------------------------------------------------------
   // Function to unroll cells dipolar field into atomic field
   //-----------------------------------------------------------------------------
   void calculate_field(const uint64_t sim_time,
                        std::vector <double>& x_spin_array, // atomic spin directions
                        std::vector <double>& y_spin_array,
                        std::vector <double>& z_spin_array,
                        std::vector <double>& m_spin_array, // atomic spin moment
                        std::vector < bool >& magnetic);    // is magnetic

   //------------------------------------------------------------------------------
   // Function to calculate energy of spin in dipole (magnetostatic) field
   //------------------------------------------------------------------------------
   double spin_magnetostatic_energy(const int atom, const double sx, const double sy, const double sz);

   //--------------------------------------------------------
   // Function to send cells field to be output in cfg file
   //--------------------------------------------------------
   int send_cells_field(std::vector<int>& cells_cell_id_array,
                        std::vector<double>& dipole_cells_field_array_x,      // B-field
                        std::vector<double>& dipole_cells_field_array_y,
                        std::vector<double>& dipole_cells_field_array_z,
                        int cells_num_local_cells
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
                   std::vector <int>& cells_num_atoms_in_cell_global, ///global  number of atoms in each cell
                   std::vector < std::vector <int> >& cells_index_atoms_array,
                   std::vector<double>& cells_volume_array,
                   std::vector<double>& cells_pos_and_mom_array,  // array to store positions and moment of cells
                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_x,
                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_y,
                   std::vector < std::vector <double> >& cells_atom_in_cell_coords_array_z,
                   std::vector<int>& atom_type_array,
                   std::vector<int>& atom_cell_id_array,
                   std::vector<double>& atom_coords_x,
                   std::vector<double>& atom_coords_y,
                   std::vector<double>& atom_coords_z,
                   std::vector<double>& x_spin_array, // atomic spin directions
                   std::vector<double>& y_spin_array,
                   std::vector<double>& z_spin_array,
                   std::vector<double>& atom_moments, // atomic magnetic moments
                   std::vector<bool>& magnetic, // bool for magnetic atoms
                   int num_atoms
   );

   //---------------------------------------------------------------------------
   // Function to process input file parameters for dipole module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   std::vector<int> get_num_atoms_in_cell_array();
   unsigned int get_num_atoms_in_cell(const int cell);
   unsigned int get_tot_num_cells();
   unsigned int get_tot_num_local_cells();
   std::vector<double> unroll_tensor(const int element, double dummy);
   std::vector<float>  unroll_tensor(const int element, float dummy);

} // end of dipole namespace

#endif //DIPOLE_H_
