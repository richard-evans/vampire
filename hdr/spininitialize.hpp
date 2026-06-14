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
   // Grain-level magnetisation modes, set via set_grain_magnetisation_mode()
   // below in response to the create:grain-magnetisation-direction keyword.
   //-----------------------------------------------------------------------------
   enum grain_magnetisation_mode_t{
      grain_mode_material    = 0, // no grain-level post-processing (default)
      grain_mode_alternating = 1  // flip spins for odd-numbered grains
   };

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
   // Function to set the grain-level magnetisation mode, called by the create
   // module when it parses the create:grain-magnetisation-direction keyword.
   // mode = 0 ("material"): no grain-level post-processing (default).
   // mode = 1 ("alternating"): the spin direction computed from the material's
   // texture is reversed (sx,sy,sz -> -sx,-sy,-sz) for every atom belonging to
   // an odd-numbered grain (grains are numbered from 0), giving alternating
   // (e.g. antiferromagnetic-like) grains.
   //---------------------------------------------------------------------------
   void set_grain_magnetisation_mode(const int mode);

   //---------------------------------------------------------------------------
   // Function to apply the grain-level magnetisation mode (see
   // set_grain_magnetisation_mode above) to the already-initialised spins of
   // every atom. This is a one-off post-processing pass over atoms::*_spin_array,
   // called once after all atoms have had their initial spin direction set by
   // initialize_spin(). Keeping this logic out of initialize_spin() avoids an
   // extra per-atom branch in the (much hotter) main initialisation loop for
   // the common case where no grain-level alternation is requested.
   //---------------------------------------------------------------------------
   void apply_grain_magnetisation_mode();

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
