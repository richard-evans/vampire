//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spintorque.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Shared variables used for the spin torque calculation
      //-----------------------------------------------------------------------------
      bool enabled=false;  // disable spin torque calculation
      bool TMRenable=false; // disable TMR calculation, GMR is set as dafault
      
      double micro_cell_size= 5*3.00; /// lateral size of spin torque microcells
      double micro_cell_thickness = 3.00; /// thickness of spin torque microcells (atomistic)

      int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      int current_direction =2; /// direction for current x->0, y->1, z->2
      //   std::vector< std::vector< micro_cell_t > > stack;
      std::vector<int> atom_st_index; // mc which atom belongs to
      std::vector<double> x_field_array; // arrays to store atomic spin torque field
      std::vector<double> y_field_array;
      std::vector<double> z_field_array;



      int num_stacks;  // total number of stacks
      int num_x_stacks; // number of stacks in x
      int num_y_stacks; // number of stack in y
      int num_microcells_per_stack; // number of microcells per stack

      int config_file_counter = 0; // spin torque config file counter

      int free_layer = 0;       /// index of free layer in magnetic tunnel junction
      int reference_layer = 0;  /// index of reference layer in magnetic tunnel junction

      double je; // = 1.0e11; // current (C/s/m^2)
      double initial_beta;
      
      double rel_angle;  // relative angle between 2 FMs for TMR calculation
      
      int ST_output_rate;
      std::vector<double> initial_m(3);

      std::vector<int> stack_index; // start of stack in microcell arrays

      std::vector<double> beta_cond; /// spin polarisation (conductivity) Beta B
      std::vector<double> beta_diff; /// spin polarisation (diffusion) Beta' Bp
      std::vector<double> sa_infinity; /// intrinsic spin accumulation
      std::vector<double> lambda_sdl; /// spin diffusion length
      std::vector<double> diffusion; /// diffusion constant Do
      std::vector<double> sd_exchange; /// electron(s)-spin(d) exchange interaction
      std::vector<double> a; // a parameter for spin accumulation
      std::vector<double> b; // b parameter for spin accumulation

      std::vector<double> coeff_ast; // adiabatic spin torque
      std::vector<double> coeff_nast; // non-adiabatic spin torque
      std::vector<double> cell_natom;


      // three-vector arrays
      std::vector<double> pos; /// stack position
      std::vector<double> m; // magnetisation
      std::vector<double> j; // spin current
      std::vector<double> sa; // spin accumulation
      std::vector<double> spin_torque; // spin torque energy (J)
      std::vector<double> ast; // adiabatic spin torque
      std::vector<double> nast; // non-adiabatic spin torque
      std::vector<double> total_ST; // non-adiabatic spin torque
      std::vector<double> magx_mat; // magnetisation of material
      std::vector<double> magy_mat;
      std::vector<double> magz_mat;
      
      
      // array of material properties
      std::vector<st::internal::mp_t> mp;

      // default material properties
      st::internal::mp_t default_properties;

   } // end of internal namespace
} // end of st namespace
