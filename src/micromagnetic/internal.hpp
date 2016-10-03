//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MICROMAGNETIC_INTERNAL_H_
#define MICROMAGNETIC_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the micromagnetic module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <vector>
#include "material.hpp"

namespace micromagnetic{

//   extern bool discretisation_micromagnetic;

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern std::vector<double> ms;
      extern std::vector<double> ku;
      extern std::vector<double> A;
      extern std::vector<double> Chi;
      extern std::vector<double> gamma;
      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      std::vector<double> calculate_ms(const int num_atoms,std::vector<double> magnetic_moment_array,std::vector<int> cell_array);
      std::vector<double> calculate_ku(const int num_atoms, const std::vector<int> cell_array, const std::vector<double> uniaxial_anisotropy_vector_x);
      std::vector<double> calculate_a(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material);
      std::vector<double> calculate_tc(int num_atoms, int num_cells, std::vector<int> cell_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material);
      std::vector<double> calculate_gamma(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array, std::vector <mp::materials_t> material);

   } // end of internal namespace

} // end of micromagnetic namespace

#endif //MICROMAGNETIC_INTERNAL_H_
