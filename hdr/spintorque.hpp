//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//   Spin torque calculation is done in 'stacks'. Each stack represents a
//   number of microcells perpendicular to the current direction
//
//                  ---------------------------
//                  | | | | | | | | | | | | | |
//                  | | | | | | | | | | | | | |
//                  | | | | | | | | | | | | | |
//                  ---------------------------
//                  | | | | | | | | | | | | | |
//           e- ->  | | | | | | | | | | | | | |
//                  | | | | | | | | | | | | | |
//                  ---------------------------
//                  | | | | | | | | | | | | | |
//                  | | | | | | | | | | | | | |
//                  | | | | | | | | | | | | | |
//                  ---------------------------
//                   ^
//                   |
//                   Starting parameters (stack[0])
//
//   Spin accumulation is calculated along each stack (in parallel). The
//   microcells in each stack are atomically thin, and definable in size.
//   The code is written so that the current direction within the sample can
//   be chosen to be either in the x,y or z direction, and the atoms are
//   assigned to cells accordingly.
//
//   The first cell in the stack contains no atoms, and is given initial
//   values for the spin transport parameters.
//
//
// System headers
#include <string>
#include <vector>

// Program headers

#ifndef SPINTORQUE_H_
#define SPINTORQUE_H_

//-----------------------------------------------------------------------------
// Namespace for variables and functions to calculate spin torque
//-----------------------------------------------------------------------------
namespace st{

   //-----------------------------------------------------------------------------
   // Variables used for the spin torque calculation
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Function to initialise spin torque calculation
   //-----------------------------------------------------------------------------
   void initialise(const double system_dimensions_x,
                  const double system_dimensions_y,
                  const double system_dimensions_z,
                  const std::vector<double>& atom_coords_x,
                  const std::vector<double>& atom_coords_y,
                  const std::vector<double>& atom_coords_z,
                  const std::vector<int>& atom_type_array,
                  const int num_local_atoms);

   //-----------------------------------------------------------------------------
   // Function to copy spin torque fields to external field array
   //-----------------------------------------------------------------------------
   void get_spin_torque_fields(std::vector<double>& x_total_external_field_array,
                               std::vector<double>& y_total_external_field_array,
                               std::vector<double>& z_total_external_field_array,
                               const int start_index,
                               const int end_index);

   //-----------------------------------------------------------------------------
   // Function for updating spin torque fields
   //-----------------------------------------------------------------------------
   void update_spin_torque_fields(const std::vector<double>& x_spin_array,
                                  const std::vector<double>& y_spin_array,
                                  const std::vector<double>& z_spin_array,
                                  const std::vector<int>& atom_type_array,
                                  const std::vector<double>& mu_s_array);

   //-----------------------------------------------------------------------------
   // Function for matching spin torque material parameters
   //-----------------------------------------------------------------------------
   bool match_material(std::string const word, std::string const value, std::string const unit, int const line, int const super_index);
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);


} // end of st namespace

#endif // SPINTORQUE_H_
