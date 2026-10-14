//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Curie Temperature program
///
/// @details Performs a temperature loop to calculate temperature dependent magnetisation
///
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.1
/// @date    09/03/2011
/// @internal
///   Created:    05/02/2011
///   Revision:   09/03/2011
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "integrate.hpp"
#include "montecarlo.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Function to calculate the temperature dependence of the anisotropy and magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature, where the constraint angles
/// are cycled. The system is initialised all spins along the constraint direction. After initialisation
/// the sytem is equilibrated for sim::equilibration timesteps before statistics are collected.
///
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    15/09/2011
///
/// @return void
///
/// @internal
///   Created:    15/09/2011
///   Revision:   --/--/----
///=====================================================================================
///
void cmc_anisotropy(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::cmc_anisotropy has been called" << std::endl;}

   // Check integrator is CMC, if not then exit disgracefully
   if(integrate::integrator != integrate::integrator_t::cmc){
      err::zexit("Program CMC-anisotropy requires Constrained Monte Carlo as the integrator. Check input file.");
   }

   // set minimum rotational angle
   // if checkpoint is loaded, then update minimum values of temperature and constrained angles
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
      zlog << zTs() << "Constrained angles theta loaded from checkpoint " << montecarlo::constraint_theta <<  std::endl;
   }
   else{
      montecarlo::constraint_theta=montecarlo::constraint_theta_min;
   }

   // perform rotational angle sweep
   while(montecarlo::constraint_theta<=montecarlo::constraint_theta_max){

      // set minimum azimuthal angle
      // if checkpoint is loaded, then update minimum values of temperature and constrained angles
      if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
         zlog << zTs() << "Constrained angle phi loaded from checkpoint: " << montecarlo::constraint_phi << std::endl;
      }
      else montecarlo::constraint_phi=montecarlo::constraint_phi_min;

      // perform azimuthal angle sweep
      while(montecarlo::constraint_phi<=montecarlo::constraint_phi_max){

         // Re-initialise spin moments for CMC
         montecarlo::CMCinit();

         // Set starting temperature
         // if checkpoint is loaded, then update minimum values of temperature
         if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag && sim::checkpoint_loaded_flag==true){
            //zlog << zTs() << "Constrained temperature loaded from checkpoint: " << sim::temperature << "K and flag = " << sim::checkpoint_loaded_flag <<std::endl;
            // In order to correctly restart the simulation, the flag that chech if the checkpoint
            // load is set to restart or load needs to be silenced
            sim::checkpoint_loaded_flag=false;
            sim::temperature += sim::delta_temperature;
         }
         else sim::temperature=sim::Tmin;

         // Perform Temperature Loop
         while(sim::temperature<=sim::Tmax){

            // Equilibrate system
            sim::integrate(sim::equilibration_time);

            // Reset mean magnetisation counters
            stats::reset();

            // Reset start time
            int start_time=sim::time;

            // Simulate system
            while(sim::time<sim::loop_time+start_time){

               // Integrate system
               sim::integrate(sim::partial_time);

               // Calculate magnetisation statistics
               stats::update();

            }

            // Output data
            vout::data();

            // Increment temperature
            sim::temperature+=sim::delta_temperature;

         } // End of temperature loop

         // Increment azimuthal angle
         montecarlo::constraint_phi+=montecarlo::constraint_phi_delta;
         montecarlo::constraint_phi_changed=true;

      } // End of azimuthal angle sweep
      if(vout::gnuplot_array_format) zmag << std::endl;

      // Increment rotational angle
      montecarlo::constraint_theta+=montecarlo::constraint_theta_delta;
      montecarlo::constraint_theta_changed=true;

   } // End of rotational angle sweep

   return;
}

}//end of namespace program
