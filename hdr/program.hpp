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

// Vampire headers
#include "program.hpp"

//==========================================================
// Namespace program
//==========================================================
namespace program
{

	//---------------------------------------------------------------------------
	// Externally visible variables
	//---------------------------------------------------------------------------
	extern int program; // program type to be run in vampire

	extern double fractional_electric_field_strength; // factor controlling strength of stt/sot and voltage

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

	// Sundry programs and diagnostics not under general release
	extern int LLB_Boltzmann();
	extern int timestep_scaling();
	extern void boltzmann_dist();
  	extern void setting_process();
  	extern void boltzmann_dist_micromagnetic_llg();

}

#endif /*PROGRAM_H_*/
