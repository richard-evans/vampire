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

#ifndef SPINTRANSPORT_INTERNAL_H_
#define SPINTRANSPORT_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the spintransport module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spintransport.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // materials class for storing exchange material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

         private:

         public:

            //----------------
            // variables
            //----------------
            double resistivity;      // spin-independent resistivity (Ohm m)
            double spin_resistivity; // spin-dependent resistivity (Ohm m)

            // constructor
            mp_t (const unsigned int max_materials = 100):
            	resistivity(1.68e-8),     // default value is for copper (Cu)
               spin_resistivity(0.0)     // default value is for copper (Cu)
            {

            }; // end of constructor

      }; // end of exchange::internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------
      extern bool enabled; // bool to enable spin transport calculation

      extern unsigned int update_rate;  // number of timesteps between updates
      extern unsigned int time_counter; // number of timesteps since last update

      extern std::vector<internal::mp_t> mp; // array of material properties

      // enumerated list of different current directions
      enum current_direction_t {px,py,pz,mx,my,mz}; // +x,+y,+z,-x,-y,-z
      extern current_direction_t current_direction; // current direction

      extern double cell_size_x; // cell size along x-direction
      extern double cell_size_y; // cell size along y-direction
      extern double cell_size_z; // cell size along z-direction

      extern unsigned int num_stacks; // number of stacks perpendicular to current direction
      extern unsigned int total_num_cells; // number of cells

      extern double voltage;                 // Applied voltage perpendicular to current direction
      extern double environment_resistivity; // Resitivity of device environment

      // array of stacks (list of cells) along current direction
      //extern std::vector < std::vector <unsigned int> > stack_array;

      // Stack data and indices
      extern std::vector <unsigned int> stack_start_index; // start of stack in 1D list of cells
      extern std::vector <unsigned int> stack_final_index; // end of stack +1 in 1D list of cells
      //extern std::vector <unsigned int> next_cell_in_stack; // list of next cell in stack to account for tunnel barrier
      extern std::vector <double> stack_resistance;        // array of stack resistances
      extern std::vector <double> stack_current;           // array of stack currents

      // arrays to store average resistance and spin resistance in each cell
      extern std::vector <double> cell_resistance;
      extern std::vector <double> cell_spin_resistance;

      // arrays to store cell properties
      extern std::vector <bool> magnetic;                    // boolean array to determine if cell is magnetic or not
      extern std::vector <double> cell_magnetization;        // 3N normalised magnetization in each cell
      extern std::vector <double> cell_isaturation;          // inverse magnetic saturation at T=0 in each cell
      extern std::vector <double> cell_position;             // 3N array of cell positions (origin)
      extern std::vector <double> cell_spin_torque_fields;   // 3N array of cell spin torque fields

      // array to store which cell each atom is in
      extern std::vector <unsigned int> atom_in_cell;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      void calculate_cell_magnetization(const unsigned int num_local_atoms,            // number of local atoms
                                        const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
                                        const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
                                        const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
                                        const std::vector<double>& atoms_m_spin_array  // moment of atoms
      );

      void calculate_magnetoresistance();

   } // end of internal namespace

} // end of spin_transport namespace

//--------------------------------------------------------------------------------
// Alias spin_transport to st namespace internal to module for code brevity
//--------------------------------------------------------------------------------
namespace st = spin_transport;

#endif //SPINTRANSPORT_INTERNAL_H_
