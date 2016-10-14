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
       std::vector <double> x_coord_array,
       std::vector <double> y_coord_array,
       std::vector <double> z_coord_array,
       double unit_cell_size_x,
       double unit_cell_size_y,
       double unit_cell_size_z,
       std::vector <double> volume_array,
       double Temperature,
       double num_atoms_in_unit_cell
   ){



      bool initialised = true;

      mm::Ax.resize(num_cells*num_cells,0.0);
      mm::Ay.resize(num_cells*num_cells,0.0);
      mm::Az.resize(num_cells*num_cells,0.0);
      mm::alpha.resize(num_cells,0.0);
      mm::chi_perp.resize(num_cells,0.0);
      mm::chi_para.resize(num_cells,0.0);
      mm::gamma.resize(num_cells,0.0);
      mm::ku.resize(num_cells,0.0);
      mm::ms.resize(num_cells,0.0);
      mm::Tc.resize(num_cells,0.0);

      mm::alpha =     mm::calculate_alpha(num_atoms, num_cells, cell_array, type_array, material);
      mm::Tc =        mm::calculate_tc(num_atoms, num_cells, cell_array,neighbour_list_array, neighbour_list_start_index, neighbour_list_end_index, type_array, material);
      mm::ku =        mm::calculate_ku(num_atoms, num_cells, cell_array, type_array, material);
      mm::gamma =     mm::calculate_gamma(num_atoms, num_cells, cell_array,type_array,material);
      mm::ms =        mm::calculate_ms(num_atoms,num_cells, cell_array, type_array,material);
      mm::chi_para =  mm::calculate_chi_para(num_cells, Temperature);
      mm::chi_perp =  mm::calculate_chi_perp(num_cells, Temperature);
      mm::Ax =         mm::calculate_ax(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array,neighbour_list_start_index,  neighbour_list_end_index,  type_array, material,unit_cell_size_x, unit_cell_size_y, unit_cell_size_z, volume_array,x_coord_array, y_coord_array, z_coord_array,num_atoms_in_unit_cell);
      mm::Ay =         mm::calculate_ay(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array,neighbour_list_start_index,  neighbour_list_end_index,  type_array, material,unit_cell_size_x, unit_cell_size_y, unit_cell_size_z, volume_array,x_coord_array, y_coord_array, z_coord_array,num_atoms_in_unit_cell);
      mm::Az =         mm::calculate_az(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array,neighbour_list_start_index,  neighbour_list_end_index,  type_array, material,unit_cell_size_x, unit_cell_size_y, unit_cell_size_z, volume_array,x_coord_array, y_coord_array, z_coord_array,num_atoms_in_unit_cell);

      std::cerr <<"\t" << mm::alpha[0] << "\t" << mm::Tc[0] << "\t" << mm::chi_para[0] << "\t" << mm::chi_perp[0] << "\t" <<mm::gamma[0] << "\t" << mm::ku[0] << "\t" << mm::ms[0] <<std::endl;

      return;

   }

} // end of micromagnetic namespace
