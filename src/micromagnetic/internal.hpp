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
#include "material.hpp"
#include "vmpi.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern int my_num_micromagnetic_cells;
      extern int my_start_index; // first cell to intergrate on local (my) cpu
      extern int my_end_index;  // last cell +1 to intergrate on local (my) cpu
      extern bool mm_correction;
      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      //vectors to store the cell parameters
      extern std::vector<double> m_e;
      extern std::vector<double> alpha_perp;
      extern std::vector<double> alpha_para;
      extern std::vector<double> A;
      extern std::vector<double> alpha;
      extern std::vector<double> one_o_chi_perp;
      extern std::vector<double> one_o_chi_para;
      extern std::vector<double> gamma;
      extern std::vector<double> ku;
      extern std::vector<double> ms;
      extern std::vector<double> Tc;

      extern std::vector <double> pinning_field_x;
      extern std::vector <double> pinning_field_y;
      extern std::vector <double> pinning_field_z;

      //stores the external fields (x,y,z)
      extern std::vector<double> ext_field;

      //stores the neighbour list for calculating A
      extern std::vector<double> macro_neighbour_list_start_index;
      extern std::vector<double> macro_neighbour_list_end_index;
      extern std::vector<double> macro_neighbour_list_array;

      extern std::vector<double> fields_neighbouring_atoms_begin;
      extern std::vector<double> fields_neighbouring_atoms_end;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      //functions to calculate the cell parameters
      std::vector<double> calculate_a(int num_atoms,
                                      int num_cells,
                                      int num_local_cells,
                                      std::vector<int> cell_array,                      //1D array storing which cell each atom is in
                                      std::vector<int> neighbour_list_array,            //1D vector listing the nearest neighbours of each atom
                                      std::vector<int> neighbour_list_start_index,      //1D vector storing the start index for each atom in the neighbour_list_array
                                      std::vector<int> neighbour_list_end_index,        //1D vector storing the end index for each atom in the neighbour_list_array
                                      const std::vector<int> type_array,                //1D array storing which material each cell is
                                      std::vector <mp::materials_t> material,           //class of material parameters for the atoms
                                      std::vector <double> volume_array,                 //1D array storing the volume of each cell
                                      std::vector <double> x_coord_array,
                                      std::vector <double> y_coord_array,
                                      std::vector <double> z_coord_array,
                                      double num_atoms_in_unit_cell,
                                      std::vector <int> local_cell_array);

      std::vector<double> calculate_alpha(int num_local_cells,int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array,
                                          std::vector <mp::materials_t> material,std::vector <int >local_cell_array);

      std::vector<double> calculate_chi_para(int num_local_cells,std::vector<int>local_cell_array,int num_cells, double T);

      std::vector<double> calculate_chi_perp(int num_local_cells,std::vector<int>local_cell_array,int num_cells, double T);

      std::vector<double> calculate_gamma(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array,
                                          std::vector <mp::materials_t> material, int num_local_cells,  std::vector<int>local_cell_array);

      std::vector<double> calculate_ku(const int num_atoms, const int num_cells, const std::vector<int> cell_array,  const std::vector<int> type_array,
                                        std::vector <mp::materials_t> material);

      std::vector<double> calculate_ms(int num_local_cells,const int num_atoms, const int num_cells,std::vector<int> cell_array, const std::vector<int> type_array,
                                        std::vector <mp::materials_t> material,std::vector <int >local_cell_array);
      std::vector<double> calculate_tc(int num_local_cells, std::vector <int> local_cell_array, int num_atoms, int num_cells, std::vector<int> cell_array, std::vector<int> neighbour_list_array,
                                      std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index,
                                      const std::vector<int> type_array, std::vector <mp::materials_t> material);

      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                                double temperature,
                                                int num_cells,
                                                int cell,
                                                std::vector<double> x_array,
                                                std::vector<double> y_array,
                                                std::vector<double> z_array);

      std::vector<double> calculate_llg_fields(std::vector <double > m,
                                                double temperature,
                                                int num_cells,
                                                int cell,
                                                std::vector<double> x_array,
                                                std::vector<double> y_array,
                                                std::vector<double> z_array);



   } // end of internal namespace

} // end of micromagnetic namespace

#endif //MICROMAGNETIC_INTERNAL_H_
