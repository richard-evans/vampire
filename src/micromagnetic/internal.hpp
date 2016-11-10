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

//#include <complex>
//#include "array3d.h"

namespace micromagnetic{

//   extern bool discretisation_micromagnetic;

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern std::vector<double> Ax;
      extern std::vector<double> Ay;
      extern std::vector<double> Az;
      extern std::vector<double> alpha;
      extern std::vector<double> chi_perp;
      extern std::vector<double> chi_para;
      extern std::vector<double> gamma;
      extern std::vector<double> ku;
      extern std::vector<double> ms;
      extern std::vector<double> Tc;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      std::vector<double> calculate_a(int num_atoms, int num_cells, int num_materials, std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material, double unit_cell_size_x, std::vector <double> volume_array, std::vector <double> x_coord_array, std::vector <double> y_coord_array, std::vector <double> z_coord_array, double num_atoms_in_unit_cell);
      std::vector<double> calculate_alpha(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array, std::vector <mp::materials_t> material);
      std::vector<double> calculate_chi_para(int num_cells, double T);
      std::vector<double> calculate_chi_perp(int num_cells, double T);
      std::vector<double> calculate_gamma(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array, std::vector <mp::materials_t> material);
      std::vector<double> calculate_ku(const int num_atoms, const int num_cells, const std::vector<int> cell_array,  const std::vector<int> type_array, std::vector <mp::materials_t> material);
      std::vector<double> calculate_ms(const int num_atoms, const int num_cells,std::vector<int> cell_array, const std::vector<int> type_array,  std::vector <mp::materials_t> material);
      std::vector<double> calculate_tc(int num_atoms, int num_cells, std::vector<int> cell_array, std::vector<int> neighbour_list_array, std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index, const std::vector<int> type_array, std::vector <mp::materials_t> material);

      void step(int num_cells, double temperature, std::vector<double> x_array,std::vector<double> y_array,std::vector<double> z_array, std::vector<double> ext_field, double dt,std::vector<double>& new_x_array,std::vector<double>& new_y_array,std::vector<double>& new_z_array);


      extern std::vector<double> ext_field;

      extern std::vector<double> x_array;
      extern std::vector<double> y_array;
      extern std::vector<double> z_array;

      extern std::vector<double> x_euler_array;
      extern std::vector<double> y_euler_array;
      extern std::vector<double> z_euler_array;

      extern std::vector<double> x_heun_array;
      extern std::vector<double> y_heun_array;
      extern std::vector<double> z_heun_array;

      extern std::vector<double> mx_store;
      extern std::vector<double> my_store;
      extern std::vector<double> mz_store;

      extern std::vector<double> mx_init;
      extern std::vector<double> my_init;
      extern std::vector<double> mz_init;


/*
      extern const double prefactor; // 1e-7/1e30

      extern Array3D<fftw_complex> Nxx; // creates the stencil complex array Nxx
      extern Array3D<fftw_complex> Nyx; // creates the stencil complex array Nyx
      extern Array3D<fftw_complex> Nzx; // creates the stencil complex array Nzx

      extern Array3D<fftw_complex> Nxy; // creates the stencil complex array Nxy
      extern Array3D<fftw_complex> Nyy; // creates the stencil complex array Nyy
      extern Array3D<fftw_complex> Nzy; // creates the stencil complex array Nzy

      extern Array3D<fftw_complex> Nxz; // creates the stencil complex array Nxz
      extern Array3D<fftw_complex> Nyz; // creates the stencil complex array Nyz
      extern Array3D<fftw_complex> Nzz; // creates the stencil complex array Nzz

      extern int num_macro_cells_x;
      extern int num_macro_cells_y;
      extern int num_macro_cells_z;
*/

   } // end of internal namespace

} // end of micromagnetic namespace

#endif //MICROMAGNETIC_INTERNAL_H_
