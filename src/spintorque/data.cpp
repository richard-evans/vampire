//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------


// C++ standard library headers
#include <vector>

// Vampire headers
#include "spintorque.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Shared variables used for the spin torque calculation
      //-----------------------------------------------------------------------------
      bool enabled=true; //false; // enable spin torque calculation

      double micro_cell_size = 5*3.54; /// lateral size of spin torque microcells
      double micro_cell_thickness = 3.54; /// thickness of spin torque microcells (atomistic)

      int current_direction; /// direction for current x->0, y->1, z->2
      //   std::vector< std::vector< micro_cell_t > > stack;
      std::vector<int> atom_st_index; // mc which atom belongs to

      int num_stacks;  // total number of stacks
      int num_x_stacks; // number of stacks in x
      int num_y_stacks; // number of stack in y
      int num_microcells_per_stack; // number of microcells per stack

      std::vector<int> stack_index; // start of stack in microcell arrays

      std::vector<double> beta_cond; /// spin polarisation (conductivity)
      std::vector<double> beta_diff; /// spin polarisation (diffusion)
      std::vector<double> sa_infinity; /// intrinsic spin accumulation
      std::vector<double> lambda_sdl; /// spin diffusion length

      // three-vector arrays
      std::vector<double> pos; /// stack position
      std::vector<double> m; // magnetisation
      std::vector<double> j; // spin current
      std::vector<double> sa; // spin accumulation
      std::vector<double> st; // spin torque
      std::vector<double> ast; // adiabatic spin torque
      std::vector<double> nast; // non-adiabatic spin torque

      // array of material properties
      std::vector<st::internal::mp_t> mp;

   } // end of internal namespace
} // end of st namespace

