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
#include <math.h>
#include "cells.hpp"
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
       double num_atoms_in_unit_cell,
       double size,
       double system_dimensions_x,
              double system_dimensions_y,
                     double system_dimensions_z
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


      mm::macro_neighbour_list_start_index.resize(num_cells,0.0);
      mm::macro_neighbour_list_end_index.resize(num_cells,0.0);



      mm::alpha =     mm::calculate_alpha(num_atoms, num_cells, cell_array, type_array, material);
      mm::Tc =        mm::calculate_tc(num_atoms, num_cells, cell_array,neighbour_list_array, neighbour_list_start_index, neighbour_list_end_index, type_array, material);
      mm::ku =        mm::calculate_ku(num_atoms, num_cells, cell_array, type_array, material);
      mm::gamma =     mm::calculate_gamma(num_atoms, num_cells, cell_array,type_array,material);
      mm::ms =        mm::calculate_ms(num_atoms,num_cells, cell_array, type_array,material);
      mm::chi_para =  mm::calculate_chi_para(num_cells, Temperature);
      mm::chi_perp =  mm::calculate_chi_perp(num_cells, Temperature);
      mm::Ax =       mm::calculate_a(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array,
         neighbour_list_start_index,  neighbour_list_end_index, type_array, material, unit_cell_size_x,
         volume_array, x_coord_array, y_coord_array, z_coord_array, num_atoms_in_unit_cell);
      //mm:: Ay =       mm::calculate_a(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array, neighbour_list_start_index,  neighbour_list_end_index, type_array, material, unit_cell_size_y,  volume_array, y_coord_array, num_atoms_in_unit_cell);
      //mm:: Az =       mm::calculate_a(num_atoms, num_cells, num_materials, cell_array, neighbour_list_array, neighbour_list_start_index,  neighbour_list_end_index, type_array, material, unit_cell_size_z,  volume_array, z_coord_array, num_atoms_in_unit_cell);

   //   for (int cell = 0; cell < num_cells; cell++) std::cerr << mm::Tc[cell] <<std::endl;

      mm::ext_field.resize(3,0.0);

      mm::x_array.resize(num_cells,0.0);
      mm::y_array.resize(num_cells,0.0);
      mm::z_array.resize(num_cells,0.0);

      mm::x_euler_array.resize(num_cells,0.0);
      mm::y_euler_array.resize(num_cells,0.0);
      mm::z_euler_array.resize(num_cells,0.0);

      mm::x_heun_array.resize(num_cells,0.0);
      mm::y_heun_array.resize(num_cells,0.0);
      mm::z_heun_array.resize(num_cells,0.0);

      mm::mx_store.resize(num_cells,0.0);
      mm::my_store.resize(num_cells,0.0);
      mm::mz_store.resize(num_cells,0.0);

      mm::mx_init.resize(num_cells,0.0);
      mm::my_init.resize(num_cells,0.0);
      mm::mz_init.resize(num_cells,0.0);

      P.resize(101);
      for (int i = 0; i < 101; ++i) P[i].resize(101);
      P1D.resize(1001,0.0);
      mean_M = 0;
      counter = 0;



//   mm::num_macro_cells_x = int(system_dimensions_x)/int(size);
//   mm::num_macro_cells_y = int(system_dimensions_y)/int(size);
//   mm::num_macro_cells_z = int(system_dimensions_z)/int(size);
/*
   std::cout << mm::num_macro_cells_x << '\t' << mm::num_macro_cells_y << '\t' << mm::num_macro_cells_z <<std::endl;

double ii,jj,kk;
   for(int i=0;i<mm::num_macro_cells_x*2;i++){
      if (i >= mm::num_macro_cells_x) ii = i - 2*mm::num_macro_cells_x;
      else ii = i;
      for(int j=0;j<mm::num_macro_cells_y*2;j++){
         if (j >= mm::num_macro_cells_y) jj = j - 2*mm::num_macro_cells_y;
         else jj = j;
         for(int k=0;k<mm::num_macro_cells_z*2;k++){
            if (k>= mm::num_macro_cells_z) kk = k - 2*mm::num_macro_cells_z;
            else kk = k;
            if((ii!=jj) && (jj != kk)){

               const double rx = ii*size; // Angstroms
               const double ry = jj*size;
               const double rz = kk*size;

               const double rij = 1.0/pow(rx*rx+ry*ry+rz*rz,0.5);

               const double ex = rx*rij;
               const double ey = ry*rij;
               const double ez = rz*rij;

               const double rij3 = rij*rij*rij; // Angstroms


               mm::Nxx(i,j,k)[0] = mm::prefactor*(3.0*ex*ex - 1.0)*rij3;
               mm::Nxy(i,j,k)[0] = mm::prefactor*(3.0*ex*ey      )*rij3;
               mm::Nxz(i,j,k)[0] = mm::prefactor*(3.0*ex*ez      )*rij3;

               mm::Nyy(i,j,k)[0] = mm::prefactor*(3.0*ey*ey - 1.0)*rij3;
               mm::Nyz(i,j,k)[0] = mm::prefactor*(3.0*ey*ez      )*rij3;
               mm::Nzz(i,j,k)[0] = mm::prefactor*(3.0*ez*ez - 1.0)*rij3;
            }
         }
      }
   }
*/
}
} // end of micromagnetic namespace
