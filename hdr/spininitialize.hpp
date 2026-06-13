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

#ifndef SPININITIALIZE_H_
#define SPININITIALIZE_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "random.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for spininitialize module
//--------------------------------------------------------------------------------
namespace spininitialize{

   //-----------------------------------------------------------------------------
   // Function to initialise spininitialize module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spininitialize module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

   //---------------------------------------------------------------------------
   // Function to get the (normalised) initial spin direction for an atom of the
   // given material at fractional coordinates (fx,fy,fz) of the system size.
   // For materials with a "random" texture, prng is used to generate the
   // random direction.
   //---------------------------------------------------------------------------
   void initialize_spin(const int material, const double fx, const double fy, const double fz,
                         double& sx, double& sy, double& sz, MTRand& prng);

   //---------------------------------------------------------------------------
   // Function to get the reference (uniform) spin direction for a material, for
   // use by programs needing a single ground-state direction per sublattice
   // (e.g. exchange stiffness calculation). Materials initialised with a
   // non-uniform texture (domain wall, skyrmion, etc.) return the default
   // direction (0,0,1).
   //---------------------------------------------------------------------------
   void get_uniform_vector(const int material, double& sx, double& sy, double& sz);

} // end of spininitialize namespace

#endif //SPININITIALIZE_H_
