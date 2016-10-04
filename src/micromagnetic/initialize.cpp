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
#include "material.hpp"

namespace mm = micromagnetic::internal;

namespace micromagnetic{

   //----------------------------------------------------------------------------
   // Function to initialize micromagnetic module
   //----------------------------------------------------------------------------
   void initialize(
       int num_cells,
       int num_atoms,
       int num_materials,
       std::vector<int> cell_array,
       std::vector<int> neighbour_list_array,
       std::vector<int> neighbour_list_start_index,
       std::vector<int> neighbour_list_end_index,
       std::vector<int> type_array,
       std::vector <mp::materials_t> material,
       class material_t
   ){



      bool initialised = true;

      mm::A.resize(num_cells*num_cells,0.0);
      mm::alpha.resize(num_cells,0.0);
      mm::chi_perp.resize(num_cells,0.0);
      mm::chi_para.resize(num_cells,0.0);
      mm::gamma.resize(num_cells,0.0);
      mm::ku.resize(num_cells,0.0);
      mm::ms.resize(num_cells,0.0);
      mm::Tc.resize(num_cells,0.0);

      std::cerr << num_cells << "\t" << num_atoms << "\t" << num_materials << std::endl;
      mm::A =         mm::calculate_a(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array,neighbour_list_start_index,  neighbour_list_end_index,  type_array, material);
      mm::alpha =     mm::calculate_alpha(num_atoms, num_cells, cell_array, type_array, material);
      mm::chi_para =  mm::calculate_chi_para(num_atoms, num_cells, cell_array, type_array, material);
      mm::chi_perp =  mm::calculate_chi_perp(num_atoms, num_cells, cell_array, type_array, material);
      mm::gamma =     mm::calculate_gamma(num_atoms, num_cells, cell_array,type_array,material);
      mm::ku =        mm::calculate_ku(num_atoms, num_cells, cell_array, type_array, material);
      mm::ms =        mm::calculate_ms(num_atoms,num_cells, cell_array, type_array,material);
      mm::Tc =        mm::calculate_tc(num_atoms, num_cells, cell_array,neighbour_list_array, neighbour_list_start_index, neighbour_list_end_index, type_array, material);


         std::cerr << num_cells << "\t" << num_atoms << "\t" << num_materials << std::endl;


      return;

   }

} // end of micromagnetic namespace
