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

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "spininitialize.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   //----------------------------------------------------------------------------
   // Function to initialize spininitialize module
   //----------------------------------------------------------------------------
   void initialize(){

      return;

   }

   //----------------------------------------------------------------------------
   // Function to get the (normalised) initial spin direction for an atom of the
   // given material at fractional coordinates (fx,fy,fz) of the system size.
   //----------------------------------------------------------------------------
   void initialize_spin(const int material, const double fx, const double fy, const double fz,
                         double& sx, double& sy, double& sz, MTRand& prng){

      // compute the spin direction from the material's texture (uniform
      // vector, random, domain wall, skyrmion, etc.)
      internal::get_spin_direction(material, fx, fy, fz, sx, sy, sz, prng);

      return;

   }

   //----------------------------------------------------------------------------
   // Function to set the grain-level magnetisation mode (see
   // spininitialize::grain_magnetisation_mode_t for the meaning of mode).
   // Called by the create module when it parses the
   // create:grain-magnetisation-direction keyword.
   //----------------------------------------------------------------------------
   void set_grain_magnetisation_mode(const int mode){

      internal::grain_magnetisation_mode = mode;

      return;

   }

   //----------------------------------------------------------------------------
   // Function to apply the grain-level magnetisation mode to the already-
   // initialised spins of every atom (see hdr/spininitialize.hpp for the full
   // rationale). Called once by create::internal::set_atom_vars after the main
   // per-atom initialisation loop has populated atoms::x/y/z_spin_array and
   // atoms::grain_array.
   //----------------------------------------------------------------------------
   void apply_grain_magnetisation_mode(){

      // default mode ("material"): spins are left exactly as computed by
      // initialize_spin(), so there is nothing to do
      if(internal::grain_magnetisation_mode == grain_mode_material) return;

      // "alternating" mode: reverse the spin direction of every atom belonging
      // to an odd-numbered grain (grains are numbered from 0, so odd grains
      // are 1, 3, 5, ...), giving neighbouring grains opposite magnetisation
      // directions (e.g. an alternating "chessboard" pattern of grains).
      for(int atom=0; atom<atoms::num_atoms; atom++){
         if(atoms::grain_array[atom] % 2 != 0){
            atoms::x_spin_array[atom] = -atoms::x_spin_array[atom];
            atoms::y_spin_array[atom] = -atoms::y_spin_array[atom];
            atoms::z_spin_array[atom] = -atoms::z_spin_array[atom];
         }
      }

      return;

   }

   //----------------------------------------------------------------------------
   // Function to get the reference (uniform) spin direction for a material
   //----------------------------------------------------------------------------
   void get_uniform_vector(const int material, double& sx, double& sy, double& sz){

      sx = 0.0;
      sy = 0.0;
      sz = 1.0;

      if(material >= 0 && material < int(internal::mp.size()) && internal::mp[material].texture == internal::uniform_vector){
         sx = internal::mp[material].initial_spin[0];
         sy = internal::mp[material].initial_spin[1];
         sz = internal::mp[material].initial_spin[2];
      }

      return;

   }

} // end of spininitialize namespace

