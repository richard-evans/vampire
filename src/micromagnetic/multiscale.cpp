//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2018. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "micromagnetic.hpp"
#include "environment.hpp"

// micromagnetic module headers
#include "internal.hpp"

//#include "random.hpp"
//#include "errors.hpp"
//#include "atoms.hpp"
#include "cells.hpp"
#include "sim.hpp"
//#include "vmpi.hpp"
//#include "vio.hpp"

namespace micromagnetic{

//------------------------------------------------------------------------------
// Serial function to compute micromagnetic and multiscale steps
//------------------------------------------------------------------------------
void multiscale_simulation_steps(const uint64_t n_steps){

   // Call environment simulation if the time is right
   if (environment::enabled && (sim::time)%environment::num_atomic_steps_env == 0){
      environment::LLB(sim::temperature, sim::H_applied,
         sim::H_vec[0], sim::H_vec[1], sim::H_vec[2], mp::dt);
   }

   // integrate for n_steps timesteps
   for(int ti=0;ti<n_steps;ti++){
      //calcaulte the field from the environment
      // if (environment::enabled && (sim::time)%environment::num_atomic_steps_env ==0) environment::LLB(sim::temperature,
      //    sim::H_applied,
      //    sim::H_vec[0],
      //    sim::H_vec[1],
      //    sim::H_vec[2],
      //    mp::dt);

      //if  there are micromagnetic cells run a micromagnetic step
      if (micromagnetic::number_of_micromagnetic_cells > 0 && (sim::time)% micromagnetic::num_atomic_steps_mm == 0) {

         //if LLG run an LLG steps
         if (micromagnetic::integrator == 0) micromagnetic::LLG(cells::local_cell_array,
            n_steps,
            cells::num_cells,
            cells::num_local_cells,
            sim::temperature,
            cells::mag_array_x,
            cells::mag_array_y,
            cells::mag_array_z,
            sim::H_vec[0],
            sim::H_vec[1],
            sim::H_vec[2],
            sim::H_applied,
            mp::dt,
            cells::volume_array);

         //if LLB run an LLB step
         else micromagnetic::LLB(cells::local_cell_array,
            n_steps,
            cells::num_cells,
            cells::num_local_cells,
            sim::temperature,
            cells::mag_array_x,
            cells::mag_array_y,
            cells::mag_array_z,
            sim::H_vec[0],
            sim::H_vec[1],
            sim::H_vec[2],
            sim::H_applied,
            mp::dt,
            cells::volume_array);
         }

         // run an atomistic step if there are atomistic atoms
         if (micromagnetic::number_of_atomistic_atoms > 0) micromagnetic::atomistic_LLG_Heun();

         //incremenet time
         sim::increment_time();

      }

      return;

   }

} // end of namespace micromagnetic
