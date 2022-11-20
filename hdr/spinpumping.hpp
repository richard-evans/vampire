//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//

#ifndef SPINPUMPING_H_
#define SPINPUMPING_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "create.hpp"
#include "spinpumping.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for spinpumping module
//--------------------------------------------------------------------------------
namespace spin_pumping{

   //-----------------------------------------------------------------------------
   // Shared variables for data output
   //-----------------------------------------------------------------------------

   // extern std::vector <double> x_s_cross_dsdt_array;        // arrays to store cross product between spin and time derivative of spins
   // extern std::vector <double> y_s_cross_dsdt_array;        // arrays to store cross product between spin and time derivative of spins
   // extern std::vector <double> z_s_cross_dsdt_array;        // arrays to store cross product between spin and time derivative of spins

   //-----------------------------------------------------------------------------
   // Function to initialise spinpumping module
   //-----------------------------------------------------------------------------
   void initialize(const double system_size_x, // maximum dimensions of system along x-direction (angstroms)
                   const double system_size_y, // maximum dimensions of system along y-direction (angstroms)
                   const double system_size_z, // maximum dimensions of system along z-direction (angstroms)
                   const int num_materials,    // number of materials
                   const uint64_t num_atoms,   // number of atoms
                   const std::vector<int>& atoms_type_array, // material types of atoms
                   const std::vector<double>& atoms_x_coord_array, // x-coordinates of atoms
                   const std::vector<double>& atoms_y_coord_array, // y-coordinates of atoms
                   const std::vector<double>& atoms_z_coord_array, // z-coordinates of atoms
                   const std::vector<double>& atoms_m_spin_array,  // moments of atoms (muB)
                   const std::vector<double>& material_damping_array, // array of material level damping constants
                   const std::vector<bool>& is_magnetic_material, // array of size num_mat to state whether material is magnetic (true) or not (false)
                   const std::vector<cs::nm_atom_t> non_magnetic_atoms_array // list of non-magnetic atoms
   );

   //-----------------------------------------------------------------------------
   // Function to update resistance, current and spin transfer torque fields
   //-----------------------------------------------------------------------------
   void update(const unsigned int num_local_atoms,            // number of local atoms
               const uint64_t time_sim,                       // simulation time
               const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
               const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
               const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
               const std::vector<double>& atoms_m_spin_array);

   void calculate_field(const unsigned int start_index,            // first atom
                        const unsigned int end_index,              // last atom
                        std::vector<double>& atoms_x_field_array,  // x-field of atoms
                        std::vector<double>& atoms_y_field_array,  // y-field of atoms
                        std::vector<double>& atoms_z_field_array); // z-field of atoms

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spinpumping module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);


} // end of spin_pumping namespace

#endif //SPINPUMPING_H_
