//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) JoelHirst 2022. All rights reserved.
//
//   Email: j.r.hirst@hallam.shu.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SPINWAVES_H_
#define SPINWAVES_H_

// C++ standard library headers
#include <string>
#include <vector> // jrh

// Vampire headers
#include "spinwaves.hpp"
#include "unitcell.hpp" // jrh
//--------------------------------------------------------------------------------
// Namespace for variables and functions for spinwaves module
//--------------------------------------------------------------------------------
namespace spinwaves{


   //-----------------------------------------------------------------------------
   // Function to initialise spinwaves module
   //-----------------------------------------------------------------------------
   extern std::vector <double> skx_r;
   extern std::vector <double> skx_i;
   extern std::vector <double> skx_r_node;
   extern std::vector <double> skx_i_node;
   extern std::vector <double> skx_r_node_transposed;
   extern std::vector <double> skx_i_node_transposed;


   // THIS IS TEMPORARY
   extern int nk_per_rank;
   extern int scatterlength;
   extern std::vector<double> skx_r_scatter;
   extern std::vector<double> skx_i_scatter;


   void fft_in_space( const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const int time );

   void fft_in_time();

   //-----------------------------------------------------------------------------
   // Function to initialise spinwaves module
   //-----------------------------------------------------------------------------
   void initialize(const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const double total_num_unit_cells_x,
                   const double total_num_unit_cells_y,
                   const double total_num_unit_cells_z,
                   const double unit_cell_size_x,
                   const double unit_cell_size_y,
                   const double unit_cell_size_z,
                   const std::vector<unitcell::atom_t>& atom,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z);

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spinwaves module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of spinwaves namespace

#endif //SPINWAVES_H_
