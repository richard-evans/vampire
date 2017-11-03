//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Adam Laverack and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "montecarlo.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for montecarlo module
//--------------------------------------------------------------------------------
namespace montecarlo{

   //-----------------------------------------------------------------------------
   // Function to initialise montecarlo module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for montecarlo module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //---------------------------------------------------------------------------
   // Function to perform one monte carlo, constrained monte carlo, or hybrid
   // cmc-mc step, respectively
   //---------------------------------------------------------------------------
   int mc_step();
   int cmc_step();
   int cmc_mc_step();

   //---------------------------------------------------------------------------
   // Provide access to CMCinit and CMCMCinit for cmc_anisotropy and
   // hybrid_cmc programs respectively
   //---------------------------------------------------------------------------
   void CMCinit();
   void CMCMCinit();

   //---------------------------------------------------------------------------
   // Function to perform monte carlo preconditioning
   //---------------------------------------------------------------------------
   void monte_carlo_preconditioning();

} // end of montecarlo namespace

#endif //MONTECARLO_H_
