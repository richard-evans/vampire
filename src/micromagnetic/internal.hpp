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
#include <iostream>
#include <fstream>  


namespace micromagnetic{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      extern int my_num_micromagnetic_cells;
      extern int my_start_index; // first cell to intergrate on local (my) cpu
      extern int my_end_index;  // last cell +1 to intergrate on local (my) cpu
      extern bool mm_correction;
      extern double res_GMR;
      extern double res_RA;
      extern std::ofstream mm_output;
      extern bool output_applied_field;
      extern bool output_m;
      extern bool output_time;
      extern bool output_temperature;
      extern bool output_resistance;
      extern std::vector <int> output_list;
      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool temperature_dependent_parameters; // flag to set temperature dependent micrmagnetic parameters

      extern std::vector < double > bias_field_x;
      extern std::vector < double > bias_field_y;
      extern std::vector < double > bias_field_z;

      extern bool bias_magnets;

      extern double bias_magnets_max_height;
      extern double bias_magnets_min_height;
      extern double shield_Ms;
      extern double bias_magnets_max_width;
      extern double bias_magnets_min_width;
      extern double bias_magnet_ms_input;
      extern double overlap_area;

      extern int bias_magnets_gap;

      extern int resistance_layer_1;
      extern int resistance_layer_2;

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
      extern std::vector<double> ku_x;
      extern std::vector<double> ku_y;
      extern std::vector<double> ku_z;
      extern std::vector<double> ms;
      extern std::vector<double> T;
      extern std::vector<double> Tc;
      extern std::vector <double> mat_vol;
      extern std::vector <double> mat_ms;
      extern std::vector <double> prefactor;

      extern std::vector <double> fmr_H;

      extern std::vector <double> pinning_field_x;
      extern std::vector <double> pinning_field_y;
      extern std::vector <double> pinning_field_z;

      extern std::vector <double> cell_material_array;

      extern double pinning_field_height;

      //stores the external fields (x,y,z)
      extern std::vector<double> ext_field;

      //stores the neighbour list for calculating A
      extern std::vector<double> macro_neighbour_list_start_index;
      extern std::vector<double> macro_neighbour_list_end_index;
      extern std::vector<double> macro_neighbour_list_array;

      extern std::vector<double> fields_neighbouring_atoms_begin;
      extern std::vector<double> fields_neighbouring_atoms_end;

      // spin transfer torque polarization vector
      extern double sttpx;
      extern double sttpy;
      extern double sttpz;

      // array to store cell level spin transfer torque parameters
      extern std::vector<double> stt_rj;
      extern std::vector<double> stt_pj;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      double calculate_resistance();
      int calculate_bias_magnets(double sx, double sy, double sz);
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

      void calculate_chi_perp(int num_mm_cells,
                             std::vector<int>& list_of_mm_cells,
                             std::vector<double>& chi_perp,
                             std::vector<double>& T,
                             std::vector<double>& Tc);

     void calculate_chi_para(int num_mm_cells,
                             std::vector<int>& list_of_mm_cells,
                             std::vector<double>& chi_para,
                             std::vector<double>& T,
                             std::vector<double>& Tc);

      std::vector<double> calculate_gamma(int num_atoms, int num_cells, std::vector<int> cell_array, const std::vector<int> type_array,
                                          std::vector <mp::materials_t> material, int num_local_cells,  std::vector<int>local_cell_array);

      std::vector<double> calculate_ku(const int num_atoms, const int num_cells, const std::vector<int> cell_array,  const std::vector<int> type_array,
                                        std::vector <mp::materials_t> material);

      // calculate spin transfer torque parameters
      void calculate_stt(const int num_atoms, const int num_cells, const std::vector<int>& cell_array, const std::vector<int>& type_array,
                         std::vector <mp::materials_t>& material, std::vector<double>& stt_rj, std::vector<double>& stt_pj);

      std::vector<double> calculate_ms(int num_local_cells,const int num_atoms, const int num_cells,std::vector<int> cell_array, const std::vector<int> type_array,
                                        std::vector <mp::materials_t> material,std::vector <int >local_cell_array);
      std::vector<double> calculate_tc(int num_local_cells, std::vector <int> local_cell_array, int num_atoms, int num_cells, std::vector<int> cell_array, std::vector<int> neighbour_list_array,
                                      std::vector<int> neighbour_list_start_index,  std::vector<int> neighbour_list_end_index,
                                      const std::vector<int> type_array, std::vector <mp::materials_t> material);

      std::vector<double> calculate_llb_fields(std::vector <double > m,
                                                double temperature,
                                                int num_cells,
                                                int cell,
                                                std::vector<double>& x_array,
                                                std::vector<double>& y_array,
                                                std::vector<double>& z_array);

      void calculate_llg_external_fields(const double temperature,
                                         const int num_cells,
                                         std::vector<double>& x_array,
                                         std::vector<double>& y_array,
                                         std::vector<double>& z_array,
                                         std::vector<double>& x_total_external_field_array,
                                         std::vector<double>& y_total_external_field_array,
                                         std::vector<double>& z_total_external_field_array);

      void calculate_llg_spin_fields(const double temperature,
                                     const int num_cells,
                                     std::vector<double>& mx_array, // Normalised magnetic moment unit vector
                                     std::vector<double>& my_array, // Normalised magnetic moment unit vector
                                     std::vector<double>& mz_array, // Normalised magnetic moment unit vector
                                     std::vector<double>& x_total_spin_field_array,   // total magnetic field
                                     std::vector<double>& y_total_spin_field_array,   // total magnetic field
                                     std::vector<double>& z_total_spin_field_array);  // total magnetic field

      void output_system_parameters(std::vector<double>& pos_and_mom_array,
                                    std::vector<double>& x_mag_array,
                                    std::vector<double>& y_mag_array,
                                    std::vector<double>& z_mag_array);

   } // end of internal namespace

} // end of micromagnetic namespace

#endif //MICROMAGNETIC_INTERNAL_H_
