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

// C++ standard library headers

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

   //----------------------------------------------------------------------------
   // Function to initialize micromagnetic module
   //----------------------------------------------------------------------------
   void initialize(
      //atom cell list
      const int num_cells,
      const int num_atoms,
      const std::vector<double> magnetic_moment_array,
      const std::vector<int> cell_array,
      const std::vector<double> uniaxial_anisotropy_vector_z,
      const std::vector<int> neighbour_list_start_index,
      const std::vector<int> neighbour_list_end_index,
      const std::vector<int> type_array,
      std::vector <mp::materials_t> material

   ){


   //   micromagnetic::ms = micromagnetic::calculate_ms(num_cells, num_atoms, magnetic_moment_array);
   //   ku = std::vector<double> calculate_ku(num_cells, num_atoms, uniaxial_anisotropy_vector_z);
   //   A = std::vector<double> calculate_A(num_cells, num_atoms);







      return;

   }

} // end of micromagnetic namespace
