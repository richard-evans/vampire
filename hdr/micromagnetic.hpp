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

// Vampire headers
#include "micromagnetic.hpp"
#include "material.hpp"
//--------------------------------------------------------------------------------
// Namespace for variables and functions for micromagnetic module
//--------------------------------------------------------------------------------
namespace micromagnetic{

   extern bool enabled; // Flag to determine if micromagnetic simulation is enabled

   // enumerated list to decide whether the simulation is micromagnetic or atomsitic
   enum discretisation_t { atomistic = 0, micromagnetics = 1, multiscale = 2 };
   extern discretisation_t discretisation_type;

   extern bool enable_resistance;
   //initialises the lists of atomstic/micromagnetic atoms for multiscale simulations
   extern std::vector < double > list_of_atomistic_atoms;
   extern std::vector < double > list_of_none_atomistic_atoms;
   extern std::vector < double > list_of_micromagnetic_cells;

   //variables to store the numbers of atomistic/ microamgnetic atoms for multiscale simulations
   extern int number_of_atomistic_atoms;
   extern int number_of_none_atomistic_atoms;
   extern int number_of_micromagnetic_cells;

   extern double MR_resistance;

   //vector to store whether cells are micromagnetic or atomistic
   extern std::vector < bool > cell_discretisation_micromagnetic;

   //set the integrator for microamgnetic simulations
   //0 - LLG
   //1 - LLB
   extern int integrator;

   //number of micromagnetic steps per atomistic step
   extern int num_atomic_steps_mm;

   //bool to enable the output of a boltzman distribution
   extern bool boltzman;

   //varibles for the Boltzman distribution
   extern double mean_M;
   extern double counter;
   extern std::vector < std::vector < double > > P;
   extern std::vector < double > P1D;

   //--------------------------------------------------------------------
   //     Function declorations
   //--------------------------------------------------------------------

   // multiscale and micromagnetic integration steps
   void multiscale_simulation_steps(const uint64_t n_steps);

   //atomsitic LLG
   int atomistic_LLG_Heun();

   //function to initialise the atomistic LLG
   int atomistic_LLGinit();

   //micromagnetic LLB
   int LLB( std::vector <int> local_cell_array,
            int num_steps,
            int num_cells,
            int num_local_cells,
            double temperature,
            std::vector<double>& x_mag_array,
            std::vector<double>& y_mag_array,
            std::vector<double>& z_mag_array,
            double Hx,
            double Hy,
            double Hz,
            double H,
            double dt,
            std::vector <double> volume_array);

    //micromagnetic LLG
    int LLG( std::vector <int> local_cell_array,
             int num_steps,
             int num_cells,
             int num_local_cells,
             double temperature,
             std::vector<double>& x_mag_array,
             std::vector<double>& y_mag_array,
             std::vector<double>& z_mag_array,
             double Hx,
             double Hy,
             double Hz,
             double H,
             double dt,
             std::vector <double> volume_array);

   //-----------------------------------------------------------------------------
   // Function to initialise micromagnetic module
   //-----------------------------------------------------------------------------
   void initialize( int num_local_cells,
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
                    std::vector <double> volume_array,
                    double Temperature,
                    double num_atoms_in_unit_cell,
                    double system_dimensions_x,
                    double system_dimensions_y,
                    double system_dimensions_z,
                    std::vector<int> local_cell_array);



   //---------------------------------------------------------------------------
   // Function to process input file parameters for micromagnetic module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of micromagnetic namespace

#endif //MICROMAGNETIC_H_
