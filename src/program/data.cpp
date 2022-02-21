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

// C++ standard library headers

// Vampire headers
#include "program.hpp"

// program module headers
#include "internal.hpp"

namespace program{

   //---------------------------------------------------------------------------
   // Externally visible variables
   //---------------------------------------------------------------------------
   int program = 0; // program type to be run in vampire
   double fractional_electric_field_strength = 1.0; // factor controlling strength of stt/sot and voltage

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside program module
      //------------------------------------------------------------------------

      bool enabled = true; // bool to enable module

      //------------------------------------------------------------------------
      // Electrial pulse program
      //------------------------------------------------------------------------
      double electrical_pulse_time      = 1.0e-9; // length of electrical pulses (1 ns default)
      double electrical_pulse_rise_time = 0.0;    // linear rise time for electrical pulse (0.0 default)
      double electrical_pulse_fall_time = 0.0;    // linear fall time for electrical pulse (0.0 default)
      int num_electrical_pulses         = 1;

      //------------------------------------------------------------------------
      // Material specific program parameters
      //------------------------------------------------------------------------
      std::vector<internal::mp_t> mp; // array of material properties

   } // end of internal namespace

} // end of program namespace
