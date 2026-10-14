//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef PROGRAM_H_
#define PROGRAM_H_

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers

//==========================================================
// Namespace program
//==========================================================
namespace program
{

   //---------------------------------------------------------------------------
   // Enumerated list of programs selected with the sim:program keyword
   //---------------------------------------------------------------------------
   // (values match the historic integer codes shown in error messages)
   enum class program_t {
      // standard programs
      benchmark = 0,
      time_series = 1,
      hysteresis_loop = 2,
      static_hysteresis_loop = 3,
      curie_temperature = 4,
      field_cool = 5,
      laser_pulse = 6,
      hamr_simulation = 7,
      cmc_anisotropy = 8,
      hybrid_cmc = 9,
      reverse_hybrid_cmc = 10,
      lagrange_multiplier = 11,
      partial_hysteresis_loop = 12,
      localised_temperature_pulse = 13,
      effective_damping = 14,
      fmr = 15,
      localised_field_cool = 16,
      electrical_pulse = 17,
      field_pulse = 18,
      laser_electrical_pulse = 19,

      // diagnostics
      diagnostic_boltzmann = 50,
      setting = 51,
      domain_wall = 52,
      exchange_stiffness = 53,
      mm_a_calculation = 54,

      // further programs not under general release
      field_sweep = 70,
      disk_tracks = 72,
      diagnostic_boltzmann_micromagnetic_llg = 73,
      spin_waves = 74,
      hamr_cool = 75
   };

   //---------------------------------------------------------------------------
   // Externally visible variables
   //---------------------------------------------------------------------------
   extern program_t program; // program type to be run in vampire

   extern double fractional_electric_field_strength; // factor controlling strength of stt/sot and voltage

   // domain wall program
   extern int domain_wall_axis;
   extern double domain_wall_position;
   extern double domain_wall_discretisation;
   extern double domain_wall_centre;
   extern double domain_wall_width;
   extern std::vector < bool > anti_PBC;
   extern std::vector < double > domain_wall_second_vector_x;
   extern std::vector < double > domain_wall_second_vector_y;
   extern std::vector < double > domain_wall_second_vector_z;

   //---------------------------------------------------------------------------
   // Function to run the program selected with the sim:program keyword
   //---------------------------------------------------------------------------
   void run();

   //---------------------------------------------------------------------------
   // Function to initialise program module
   //---------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for program module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   // program functions
   extern int bmark();
   extern void time_series();
   extern int hysteresis();
   extern int static_hysteresis();
   extern void partial_hysteresis_loop();
   extern int curie_temperature();
   extern void field_cool();
   extern void local_field_cool();
   extern void temperature_pulse();
   extern void hamr();
   extern void cmc_anisotropy();
   extern void hybrid_cmc();
   extern void reverse_hybrid_cmc();
   extern void lagrange_multiplier();
   extern void localised_temperature_pulse();
   extern void effective_damping();
   extern void fmr();
   extern void tracks();
   extern void field_sweep();
   extern void fmr();
   extern void domain_wall();
   extern void mm_A_calculation();
   extern void exchange_stiffness();
   extern void electrical_pulse();
   extern void laser_electrical_pulse();
   extern void spin_waves(); // JRH
   extern void field_pulse();

   // Sundry programs and diagnostics not under general release
   extern int LLB_Boltzmann();
   extern int timestep_scaling();
   extern void boltzmann_dist();
   extern void setting_process();
   extern void boltzmann_dist_micromagnetic_llg();
   extern void hamr_cool();

}

#endif /*PROGRAM_H_*/
