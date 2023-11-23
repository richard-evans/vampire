//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SLD_H_
#define SLD_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "sld.hpp"
#include <vector>


//--------------------------------------------------------------------------------
// Namespace for variables and functions for sld module
//--------------------------------------------------------------------------------
namespace sld{

   //-----------------------------------------------------------------------------
   // Function to initialise sld module
   //-----------------------------------------------------------------------------
   void initialize();
   void tests();


   //---------------------------------------------------------------------------
   // Function to process input file parameters for sld module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //---------------------------------------------------------------------------
   // Function to calculate forces in the spin-lattice module
   //---------------------------------------------------------------------------
   void compute_forces(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& coord_array_x0, // coord vectors for atoms
               const std::vector<double>& coord_array_y0,
               const std::vector<double>& coord_array_z0,
               const std::vector<double>& coord_array_x, // coord vectors for atoms
               const std::vector<double>& coord_array_y,
               const std::vector<double>& coord_array_z,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& potential_eng);
   //
   void compute_fields(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& coord_array_x, // coord vectors for atoms
               const std::vector<double>& coord_array_y,
               const std::vector<double>& coord_array_z,
               const std::vector<double>& spin_array_x, // coord vectors for atoms
               const std::vector<double>& spin_array_y,
               const std::vector<double>& spin_array_z,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& fields_array_x, //  vectors for forces
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z);

   void compute_forces_fields(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& coord_array_x0, // coord vectors for atoms
               const std::vector<double>& coord_array_y0,
               const std::vector<double>& coord_array_z0,
               std::vector<double>& coord_array_x, // coord vectors for atoms
               std::vector<double>& coord_array_y,
               std::vector<double>& coord_array_z,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z);

   extern int suzuki_trotter();

   extern void stats_sld();

   extern std::vector<double> forces_array_x;
   extern double var_test;
   double PBC_wrap ( double dx, double L, bool bounds);

   double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               std::vector<double>& x_spin_array, // coord vectors for atoms
               std::vector<double>& y_spin_array,
               std::vector<double>& z_spin_array,
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z,
               std::vector<double>& mu_s_array);

   double compute_lattice_temperature(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               std::vector<double>& velo_array_x, // coord vectors for atoms
               std::vector<double>& velo_array_y,
               std::vector<double>& velo_array_z);

//
   double compute_potential_energy(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array);

   double compute_kinetic_energy(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               std::vector<double>& velo_array_x, // coord vectors for atoms
               std::vector<double>& velo_array_y,
               std::vector<double>& velo_array_z);

//
    double compute_effective_J(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            std::vector<double>& sum_J);
//  //
    double compute_effective_C(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            std::vector<double>& sum_C);
//
   double compute_exchange_energy(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index);
//
   double compute_coupling_energy(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index);

//
   extern double lattice_temperature;
   extern double spin_temperature;

/*
   extern double potential_energy;
   extern double kinetic_energy;
   extern double sld_exchange_energy;
   extern double sld_coupling_energy;
   extern double sld_total_energy;
   extern double sld_total_spin_energy;*/  
   
   
   extern double J_eff;
   extern double C_eff;
   
   void suzuki_trotter_parallel_init(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                         double min_dim[3], double max_dim[3]);
   extern bool suzuki_trotter_parallel_initialized;
   void suzuki_trotter_step_parallel(std::vector<double> &x_spin_array, std::vector<double> &y_spin_array, std::vector<double> &z_spin_array, std::vector<int> &type_array);



} // end of sld namespace

#endif //SLD_H_
