//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack Collings and Sam Westmoreland 2016. All rights reserved.
//
//   Email: sw766@york.ac.uk
//          jbc525@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ANISOTROPY_H_
#define ANISOTROPY_H_

// C++ standard library headers
//#include <vector>
//#include <string>

// Vampire headers
#include "anisotropy.hpp"
#include "create.hpp"

//--------------------------------------------------------------------------------
// namespace for variables and functions for anisotropy module
//--------------------------------------------------------------------------------
namespace anisotropy
{
   //-----------------------------------------------------------------------------
   // Function to check orthogonality, set orthonormality and generate last
   // last anisotropy vector
   //-----------------------------------------------------------------------------
   void set_anisotropy_vectors(  std::vector<double>& v1, 
                                 std::vector<double>& v2,
                                 std::vector<double>& v3,
                                 int mat);

   //-----------------------------------------------------------------------------
   // function to initialize anisotropy module
   //-----------------------------------------------------------------------------
   void initialize(const unsigned int num_atoms,
                   std::vector<int>& atom_type_array,
                   std::vector<double>& mu_s_array);

   //-----------------------------------------------------------------------------
   // function to calculate anisotropy fields
   //-----------------------------------------------------------------------------
   void fields(std::vector<double>& spin_array_x,
               std::vector<double>& spin_array_y,
               std::vector<double>& spin_array_z,
               std::vector<int>&    type_array,
               std::vector<double>& field_array_x,
               std::vector<double>& field_array_y,
               std::vector<double>& field_array_z,
               const int start_index,
               const int end_index,
               const double temperature);

   //-----------------------------------------------------------------------------
   // function to calculate anisotropy energy for a single spin
   //-----------------------------------------------------------------------------
   double single_spin_energy(const int atom, const int imaterial, const double sx, const double sy, const double sz, const double temperature);

   //-----------------------------------------------------------------------------
   // functions to get anisotropy constants
   //-----------------------------------------------------------------------------
   double get_anisotropy_constant(const int material);
   double get_ku2(const int material);
   double get_k2r1(const int material);
   double get_k2r1_odd(const int material);
   double get_k2r2(const int material);
   double get_k2r2_odd(const int material);
   double get_ku4(const int material);
   double get_k4r1(const int material);
   double get_k4r1_odd(const int material);
   double get_k4r2(const int material);
   double get_k4r2_odd(const int material);
   double get_k4r3(const int material);
   double get_k4r3_odd(const int material);
   double get_k4r4(const int material);
   double get_k4r4_odd(const int material);
   double get_ku6(const int material);
   double get_k6r1_odd(const int material);
   double get_k6r2(const int material);
   double get_k6r2_odd(const int material);
   double get_k6r3(const int material);
   double get_k6r3_odd(const int material);
   double get_k6r4(const int material);
   double get_k6r4_odd(const int material);
   double get_k6r5_odd(const int material);
   double get_k6r5(const int material);
   double get_k6r6(const int material);
   double get_k6r6_odd(const int material);
   double get_kc4(const int material);
   double get_kc6(const int material);
   std::vector<double> get_ku_vector(const int material);
   std::vector<double> get_kr_vector(const int material);
   std::vector<double> get_kl_vector(const int material);

   //-----------------------------------------------------------------------------
   // function to identify surface atoms
   //-----------------------------------------------------------------------------
   void identify_surface_atoms(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <neighbours::neighbour_t> > & cneighbourlist);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for anisotropy module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, const int max_materials);

} // end of anisotropy namespace

#endif //ANISOTROPY_H_
