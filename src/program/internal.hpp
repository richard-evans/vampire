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

#ifndef PROGRAM_INTERNAL_H_
#define PROGRAM_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the program module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <string>
#include <vector>

// Vampire headers
#include "program.hpp"

// program module headers
#include "internal.hpp"

namespace program {

namespace internal {

//-------------------------------------------------------------------------
// Internal data type definitions
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// internal materials class for storing material parameters
//-----------------------------------------------------------------------------
class mp_t {

private:
public:
  //------------------------------
  // material parameter variables
  //------------------------------
  double test;

  // constructor
  mp_t(const unsigned int max_materials = 100)
      : test(0.0) // constructor initialisation of test variable
  {
    // constructor body for initialising more complex data/arrays
  }; // end of constructor

}; // end of internal::mp class

//-------------------------------------------------------------------------
// Internal shared variables
//-------------------------------------------------------------------------
extern bool enabled; // bool to enable module

//------------------------------------------------------------------------
// Electrial pulse program
//------------------------------------------------------------------------
extern double
    electrical_pulse_time; // length of electrical pulses (1 ns default)
extern double electrical_pulse_rise_time; // linear rise time for electrical
                                          // pulse (0.0 default)
extern double electrical_pulse_fall_time; // linear fall time for electrical
                                          // pulse (0.0 default)
extern double electrical_strength;        // pulse strength

extern double electrical_duration;

extern double electrical_begin_time;

extern double electrical_consistent_time;

extern std::string electrical_pulse_shape; // pulse shape (square and gaussian)

extern int num_electrical_pulses;

//------------------------------------------------------------------------
// Field pulse program
//------------------------------------------------------------------------
extern double field_pulse_time; // length of field pulses (1 ns default)

//------------------------------------------------------------------------
// Exchange stiffness program
//------------------------------------------------------------------------
extern double exchange_stiffness_max_constraint_angle;   // degrees
extern double exchange_stiffness_delta_constraint_angle; // degrees

//------------------------------------------------------------------------
// Material level parameters
//------------------------------------------------------------------------
extern std::vector<internal::mp_t> mp; // array of material properties

//-------------------------------------------------------------------------
// Internal function declarations
//-------------------------------------------------------------------------

} // namespace internal

} // namespace program

#endif // PROGRAM_INTERNAL_H_
