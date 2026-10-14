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

   // LaGrange multiplier program
   double lagrange_lambda_x=0.0;
   double lagrange_lambda_y=0.0;
   double lagrange_lambda_z=0.0;
   double lagrange_m=1.0;
   double lagrange_N=1000.0;
   bool   lagrange_multiplier_enabled=false;

   // disk tracks program (track and LFA parameters)
   std::vector < double > track_field_x;
   std::vector < double > track_field_y;
   std::vector < double > track_field_z;

   // laser pulse programs (pump and two-temperature model parameters)
   double cooling_time=100.0e-12; ///seconds
   int cooling_function_flag=0; /// 0 = exp, 1 = gaussian
   pump_functions_t pump_function=two_temperature;
   double pump_power=20.0; // mJ/cm^2;
   double pump_time=50.0e-15;
   double double_pump_power=20.0; // mJ/cm^2;
   double double_pump_Tmax=500.0;
   double double_pump_time=50.0e-15;
   double double_pump_delay=10.0e-12;
   double HeatSinkCouplingConstant=0.0; ///1.1e12 ~ sensible value
   double TTCe = 222.0; ///electron specific heat (gamma)
   double TTCl = 2.3E06; ///phonon specific heat
   double TTG = 6.6E17 ;///electron coupling constant
   double TTTe = 0.0; /// electron temperature
   double TTTp = 0.0; /// phonon temperature

   double track_Ms = 0.1;
   double track_bit_width = 800;
   double track_bit_size = 1000;
   double track_bit_depth = 600;
   double cross_track_velocity = 0.0;
   double down_track_velocity = 0.0;
   double track_pos_x;
   double track_pos_z;
   double Ms;


   double initial_down_track_position = 0.0;
   double initial_cross_track_position = 0.0;

   int track_num_bits_per_track = 1;
   int track_num_tracks = 1;

   double track_bit_gap = 10.0;
    double track_track_gap = 10.0;

   // distance of tracks from read head
   double track_fly_height = 100.0; // Angstroms

   double LFA_scan_field_step = 0.01;
   bool LFA = false;

   bool track_ms_file = false;

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
