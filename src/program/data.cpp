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
   program_t program = program_t::benchmark; // program type to be run in vampire
   double fractional_electric_field_strength = 1.0; // factor controlling strength of stt/sot and voltage

   // domain wall program
   int domain_wall_axis = 0;
   double domain_wall_position = 0.25;
   double domain_wall_discretisation = 10;
   double domain_wall_centre = 0;
   double domain_wall_width = 10.0;
   std::vector <bool > anti_PBC(3,false);
   std::vector < double > domain_wall_second_vector_x(100,0);
   std::vector < double > domain_wall_second_vector_y(100,0);
   std::vector < double > domain_wall_second_vector_z(100,1.0);

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
      double electrical_pulse_delay     = 0.0;    // delay from laser pulse centre to electrical pulse centre (0.0 default)

      //------------------------------------------------------------------------
      // Field pulse program
      //------------------------------------------------------------------------
      double field_pulse_time      = 1.0e-9; // length of field pulses (1 ns default)

      //------------------------------------------------------------------------
      // Exchange stiffness program
      //------------------------------------------------------------------------
      double exchange_stiffness_max_constraint_angle   = 180.01; // degrees
      double exchange_stiffness_delta_constraint_angle =  5; // 22.5 degrees

      //------------------------------------------------------------------------
      // Effective damping program
      //------------------------------------------------------------------------
      double effective_damping_rotation_angle = 30.0; // initial rotation of magnetisation from z-axis (degrees)

      //------------------------------------------------------------------------
      // Material specific program parameters
      //------------------------------------------------------------------------
      std::vector<internal::mp_t> mp; // array of material properties

   } // end of internal namespace

} // end of program namespace
