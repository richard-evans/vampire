//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "sim.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{
   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------

   int num_monte_carlo_preconditioning_steps(0);

   uint64_t time         = 0; // time step counter
   uint64_t total_time   = 10000; // total time steps (non-loop code)
   uint64_t loop_time    = 10000; // loop time steps (hysteresis/temperature loops)
   uint64_t partial_time = 1000; // same as time-step-increment
   uint64_t equilibration_time = 0; // equilibration time steps

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
      //----------------------------------------------------------------------------
      // Shared variables used within sim module
      //---------------------------------------------------------------------------
      std::vector<sim::internal::mp_t> mp; // array of material properties
      std::vector<double> slonczewski_aj; // array of adiabatic spin torques
      std::vector<double> slonczewski_bj; // array of non-adiabatic spin torques
      std::vector<double> slonczewski_spin_polarization_unit_vector(3,0.0); // spin polarization direction
      std::vector<double> SOT_DL; // array of adiabatic spin torques
      std::vector<double> SOT_FL; // array of non-adiabatic spin torques
      std::vector<double> SOT_spin_polarization_unit_vector(3,0.0); // spin polarization direction


   } // end of internal namespace

} // end of sim namespace
