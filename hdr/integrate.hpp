//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

// C++ standard library headers
#include <cstdint>
#include <string>

//--------------------------------------------------------------------------------
// Namespace for variables and functions for integrate module
//--------------------------------------------------------------------------------
namespace integrate{

   //-----------------------------------------------------------------------------
   // Enumerated list of integrators selected with the sim:integrator keyword
   //-----------------------------------------------------------------------------
   // (values match the historic integer codes shown in error messages)
   enum class integrator_t {
      // LLG spin dynamics
      llg_heun = 0,
      llg_midpoint = 2,
      llg_quantum = 5,
      llg_rkmk2 = 10,
      llg_rkmk4 = 11,
      llg_dm = 12,

      // longitudinal spin fluctuations
      lsf = 6,
      lsf_rk4 = 8,

      // Monte Carlo (montecarlo module)
      monte_carlo = 1,
      cmc = 3,
      hybrid_cmc = 4,
      lsf_mc = 7,

      // spin-lattice dynamics (sld module)
      suzuki_trotter = 9
   };

   //-----------------------------------------------------------------------------
   // Externally visible variables
   //-----------------------------------------------------------------------------
   extern integrator_t integrator; // variable to specify integrator

   //-----------------------------------------------------------------------------
   // Function to integrate the system for n_steps time steps (serial or MPI)
   //-----------------------------------------------------------------------------
   int integrate(uint64_t n_steps);

   //-----------------------------------------------------------------------------
   // Legacy Landau-Lifshitz-Bloch integrator (used by LLB-Boltzmann program)
   //-----------------------------------------------------------------------------
   int LLB(const int num_steps);

   //-----------------------------------------------------------------------------
   // Function to initialise integrate module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for integrate module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to get LSF Landau coefficients (J) for a material
   //---------------------------------------------------------------------------
   void get_lsf_coefficients(const int material, double& A, double& B, double& C);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of integrate namespace

#endif //INTEGRATE_H_
