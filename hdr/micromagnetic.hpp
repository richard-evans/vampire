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

#ifndef MICROMAGNETIC_H_
#define MICROMAGNETIC_H_

// C++ standard library headers
#include <string>
#include <vector>
// Vampire headers
#include "micromagnetic.hpp"
#include "material.hpp"
//--------------------------------------------------------------------------------
// Namespace for variables and functions for micromagnetic module
//--------------------------------------------------------------------------------
namespace micromagnetic{

   extern bool discretisation_micromagnetic;
   extern bool initialised;


   void initialize(int num_cells,int num_atoms,int num_materials, std::vector<int> cell_array,std::vector<int> neighbour_list_array,std::vector<int> neighbour_list_start_index, std::vector<int> neighbour_list_end_index,std::vector<int> type_array,std::vector <mp::materials_t> material,std::vector <double> x_coord_array,std::vector <double> y_coord_array, std::vector <double> z_coord_array, double unit_cell_size_x, double unit_cell_size_y, double unit_cell_size_z, std::vector <double> volume_array, double Temperature, double num_atoms_in_unit_cell, double size, double system_dimensions_x,double system_dimensions_y,double system_dimensions_z );
   extern int LLB(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H, double dt, std::vector <double> volume_array, int N);
   extern bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);
   extern bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   extern std::vector < std::vector<int > > P;
   extern std::vector < int > P1D;
   extern double mean_M;
   extern int counter;


} // end of micromagnetic namespace

#endif //MICROMAGNETIC_H_
