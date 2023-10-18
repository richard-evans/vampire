//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cmath>
#include <vector>

// Vampire Header files
#include "random.hpp"
#include "sim.hpp"
#include "atoms.hpp"
#include "material.hpp"

// Internal header
#include "internal.hpp"

namespace montecarlo{
//------------------------------------------------------------------------------
// Integrates a Monte Carlo step
//------------------------------------------------------------------------------
void mc_step(std::vector<double> &x_spin_array,
             std::vector<double> &y_spin_array,
             std::vector<double> &z_spin_array,
             int num_atoms,
             std::vector<int> &type_array){

      // calculate number of steps to calculate
      const int nmoves = num_atoms;

      // Temporaries
      int atom=0;
      double Eold=0.0;
      double Enew=0.0;
      double DE=0.0;

      // Material dependent temperature rescaling
      std::vector<double> rescaled_material_kBTBohr(internal::num_materials);
      std::vector<double> sigma_array(internal::num_materials); // range for tuned gaussian random move
      for(int m=0; m<internal::num_materials; ++m){
         double alpha = internal::temperature_rescaling_alpha[m];
         double Tc = internal::temperature_rescaling_Tc[m];
         double rescaled_temperature = sim::temperature < Tc ? Tc*pow(sim::temperature/Tc,alpha) : sim::temperature;
         rescaled_material_kBTBohr[m] = 9.27400915e-24/(rescaled_temperature*1.3806503e-23);
         sigma_array[m] = rescaled_temperature < 1.0 ? 0.02 : pow(1.0/rescaled_material_kBTBohr[m],0.2)*0.08;
      }

      double statistics_moves = 0.0;
      double statistics_reject = 0.0;

      // loop over natoms to form a single Monte Carlo step
      for(int i=0;i<nmoves; i++){

         // add one to number of moves counter
         statistics_moves+=1.0;

         // pick random atom
         atom = int(nmoves*mtrandom::grnd());

         // get material id
         const int imaterial=type_array[atom];

         // Calculate range for move
         internal::delta_angle=sigma_array[imaterial];

         // Save old spin position
         internal::Sold[0] = x_spin_array[atom];
         internal::Sold[1] = y_spin_array[atom];
         internal::Sold[2] = z_spin_array[atom];

         // Make Monte Carlo move
         internal::mc_move(internal::Sold, internal::Snew);

         // Calculate current energy
         Eold = sim::calculate_spin_energy(atom);

         // Copy new spin position
         x_spin_array[atom] = internal::Snew[0];
         y_spin_array[atom] = internal::Snew[1];
         z_spin_array[atom] = internal::Snew[2];

         // Calculate new energy
         Enew = sim::calculate_spin_energy(atom);

         // Calculate difference in Joules/mu_B
         DE = (Enew-Eold)*internal::mu_s_SI[imaterial]*1.07828231e23; //1/9.27400915e-24

         // Check for lower energy state and accept unconditionally
         if(DE<0) continue;
         // Otherwise evaluate probability for move
         else{
            if(exp(-DE*rescaled_material_kBTBohr[imaterial]) >= mtrandom::grnd()) continue;
            // If rejected reset spin coordinates and continue
            else{
               x_spin_array[atom] = internal::Sold[0];
               y_spin_array[atom] = internal::Sold[1];
               z_spin_array[atom] = internal::Sold[2];
               // add one to rejection counter
               statistics_reject += 1.0;
               continue;
            }
         }
      }

      // calculate new adaptive step sigma angle
      if(montecarlo::algorithm == montecarlo::adaptive){
         const double last_rejection_rate = statistics_reject / statistics_moves;
         const double factor = 0.5 / last_rejection_rate;
         montecarlo::internal::adaptive_sigma *= factor;
         // check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
      	if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5) montecarlo::internal::adaptive_sigma = 60.0;
      }

      // Save statistics to sim namespace variable
      sim::mc_statistics_moves += statistics_moves;
      sim::mc_statistics_reject += statistics_reject;

      return;

   }

   } // End of namespace montecarlo
