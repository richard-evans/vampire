//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spintransport.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   double total_resistance = 0.0;
   double total_current = 0.0;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside spintransport module
      //------------------------------------------------------------------------
      bool enabled = false; // bool to enable spin transport calculation

      unsigned int update_rate  = 100;  // number of timesteps between updates
      unsigned int time_counter = 100;  // number of timesteps since last update (initially set to same as update rate to ensure calculation at start)

      std::vector<internal::mp_t> mp; // array of material properties

      // enumerated list of different current directions
      current_direction_t current_direction = st::internal::pz; // current direction (default along +z direction)

      double cell_size_x = 10.0; // cell size along x-direction (1 nm default size)
      double cell_size_y = 10.0; // cell size along y-direction (1 nm default size)
      double cell_size_z = 10.0; // cell size along z-direction (1 nm default size)

      unsigned int num_stacks = 0; // number of stacks perpendicular to current direction
      unsigned int total_num_cells  = 0; // number of cells

      double voltage = 0.0;                    // Applied voltage perpendicular to current direction
      double environment_resistivity = 1.0e15; // Resitivity of device environment (assumed to be SiO2)

      // Stack data and indices
      std::vector <unsigned int> stack_start_index; // start of stack in 1D list of cells
      std::vector <unsigned int> stack_final_index; // end of stack +1 in 1D list of cells
      std::vector <double> stack_resistance;        // array of stack resistances
      std::vector <double> stack_current;           // array of stack currents

      unsigned int first_stack = 0; // first stack on local processor
      unsigned int last_stack = 0;  // last stack on local processor

      // arrays to store average resistance and spin resistance in each cell
      std::vector <double> cell_resistance;
      std::vector <double> cell_spin_resistance;

      // arrays to store cell properties
      std::vector <bool> magnetic;                    // boolean array to determine if cell is magnetic or not
      std::vector <double> cell_alpha;                // cell magnetization (average of constituent atoms)
      std::vector <double> cell_magnetization;        // 3N normalised magnetization in each cell
      std::vector <double> cell_isaturation;          // inverse magnetic saturation at T=0 in each cell
      std::vector <double> cell_position;             // 3N array of cell positions (origin)
      std::vector <double> cell_spin_torque_fields;   // 3N array of cell spin torque fields
      std::vector <double> cell_relaxation_torque_rj; // cell specific prefactors for spin-torque relaxation bj
      std::vector <double> cell_precession_torque_pj; // cell specific prefactors for spin-torque precession aj

      // array to store which cell each atom is in
      std::vector <unsigned int> atom_in_cell;

   } // end of internal namespace

} // end of spin_transport namespace
