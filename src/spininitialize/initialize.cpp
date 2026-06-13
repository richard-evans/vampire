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

      internal::get_spin_direction(material, fx, fy, fz, sx, sy, sz, prng);

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

