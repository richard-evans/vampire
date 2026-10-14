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
#include "material.hpp"
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
///   Created:    02/11/2011
///   Revision:   --/--/----
///=====================================================================================
///
void hybrid_cmc(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::hybrid_cmc has been called" << std::endl;}

   // Check integrator is CMC, if not then exit disgracefully
   if(integrate::integrator != integrate::integrator_t::hybrid_cmc){
      err::fatal("Error! cmc-anisotropy program requires Hybrid Constrained Monte Carlo as the integrator. Exiting.");
   }

   // resize cmc array to include correct number of materials
   montecarlo::cmc::cmc_mat.resize(mp::num_materials);

   // loop over all materials
   for (int mat=0;mat<mp::num_materials;mat++){

   // Set active material (used for calculating rotational update of spin directions)
   montecarlo::cmc::active_material=mat;

   // Output informative message to log file
   zlog << zTs() << "Starting Hybrid CMC loop for material " << mat << std::endl;

   // set minimum rotational angle
   montecarlo::cmc::cmc_mat[mat].constraint_theta=montecarlo::cmc::cmc_mat[mat].constraint_theta_min;

   // perform rotational angle sweep
   while(montecarlo::cmc::cmc_mat[mat].constraint_theta<=montecarlo::cmc::cmc_mat[mat].constraint_theta_max){

      // set minimum azimuthal angle
      montecarlo::cmc::cmc_mat[mat].constraint_phi=montecarlo::cmc::cmc_mat[mat].constraint_phi_min;

      // perform azimuthal angle sweep
      while(montecarlo::cmc::cmc_mat[mat].constraint_phi<=montecarlo::cmc::cmc_mat[mat].constraint_phi_max){

         // Re-initialise spin moments for CMC
         montecarlo::CMCMCinit();

         // Set starting temperature
         sim::temperature=sim::Tmin;

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
         if(montecarlo::cmc::cmc_mat[mat].constraint_phi+montecarlo::cmc::cmc_mat[mat].constraint_phi_delta>montecarlo::cmc::cmc_mat[mat].constraint_phi_max) break;
         montecarlo::cmc::cmc_mat[mat].constraint_phi+=montecarlo::cmc::cmc_mat[mat].constraint_phi_delta;
         montecarlo::constraint_phi_changed=true;
      } // End of azimuthal angle sweep
      if(vout::gnuplot_array_format) zmag << std::endl;

      // Increment rotational angle
      if(montecarlo::cmc::cmc_mat[mat].constraint_theta+montecarlo::cmc::cmc_mat[mat].constraint_theta_delta>montecarlo::cmc::cmc_mat[mat].constraint_theta_max) break;
      montecarlo::cmc::cmc_mat[mat].constraint_theta+=montecarlo::cmc::cmc_mat[mat].constraint_theta_delta;
      montecarlo::constraint_theta_changed=true;
   } // End of rotational angle sweep

   } // end of material loop
   return;
}

/// @brief Function to calculate the temperature dependence of the anisotropy and magnetisation
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature, where the constraint angles
/// are cycled. The system is initialised all spins along the constraint direction. After initialisation
/// the sytem is equilibrated for sim::equilibration timesteps before statistics are collected. In this
/// version the order of constrained angles is changed.
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
///   Created:    02/11/2011
///   Revision:   --/--/----
///=====================================================================================
///
void reverse_hybrid_cmc(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "program::reversed_hybrid_cmc has been called" << std::endl;}

   // Check integrator is CMC, if not then exit disgracefully
   if(integrate::integrator != integrate::integrator_t::hybrid_cmc){
     err::fatal("Error! cmc-anisotropy program requires Hybrid Constrained Monte Carlo as the integrator. Exiting.");
   }

   // resize cmc array to include correct number of materials
   montecarlo::cmc::cmc_mat.resize(mp::num_materials);

   // loop over all materials
   for (int mat=0;mat<mp::num_materials;mat++){

   // Set active material (used for calculating rotational update of spin directions)
   montecarlo::cmc::active_material=mat;

   // Output informative message to log file
   zlog << zTs() << "Starting Hybrid CMC loop for material " << mat << std::endl;

   // set minimum azimuthal angle
   montecarlo::cmc::cmc_mat[mat].constraint_phi=montecarlo::cmc::cmc_mat[mat].constraint_phi_min;

   // perform azimuthal angle sweep
   while(montecarlo::cmc::cmc_mat[mat].constraint_phi<=montecarlo::cmc::cmc_mat[mat].constraint_phi_max){

      // set minimum rotational angle
      montecarlo::cmc::cmc_mat[mat].constraint_theta=montecarlo::cmc::cmc_mat[mat].constraint_theta_min;

      // perform rotational angle sweep
      while(montecarlo::cmc::cmc_mat[mat].constraint_theta<=montecarlo::cmc::cmc_mat[mat].constraint_theta_max){

         // Re-initialise spin moments for CMC
         montecarlo::CMCMCinit();

         // Set starting temperature
         sim::temperature=sim::Tmin;

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

         // Increment rotational angle
         if(montecarlo::cmc::cmc_mat[mat].constraint_theta+montecarlo::cmc::cmc_mat[mat].constraint_theta_delta>montecarlo::cmc::cmc_mat[mat].constraint_theta_max) break;
         montecarlo::cmc::cmc_mat[mat].constraint_theta+=montecarlo::cmc::cmc_mat[mat].constraint_theta_delta;
         montecarlo::constraint_theta_changed=true;
      } // End of rotational angle sweep
      if(vout::gnuplot_array_format) zmag << std::endl;

      // Increment azimuthal angle
      if(montecarlo::cmc::cmc_mat[mat].constraint_phi+montecarlo::cmc::cmc_mat[mat].constraint_phi_delta>montecarlo::cmc::cmc_mat[mat].constraint_phi_max) break;
      montecarlo::cmc::cmc_mat[mat].constraint_phi+=montecarlo::cmc::cmc_mat[mat].constraint_phi_delta;
      montecarlo::constraint_phi_changed=true;
   } // End of azimuthal angle sweep

   } // end of material loop
   return;
}

}//end of namespace program
