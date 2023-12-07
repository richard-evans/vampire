// Standard libraries
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>

// Vampire header files
#include "../simulate/internal.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "exchange.hpp"
#include "internal.hpp"
#include "lsf_mc.hpp"
#include "material.hpp"
#include "montecarlo.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmath.hpp"

namespace montecarlo{

   std::vector<double> mod_S;
   bool mc_set = false;

   void mcinit(){

      // Initialize spin length to 1.0
      mod_S.resize(atoms::num_atoms, 1.0);

      mc_set = true;
   }

   /// Spin move
   /// Move spin within cone near old position
   void mc_angle_transverse(const std::vector<double> &old_spin, std::vector<double> &new_spin, const double angle, const int atom){

      double factorx = mtrandom::grnd();
      double factory = mtrandom::grnd();
      double factorz = mtrandom::grnd();

      // Function to ensure random numbers have a 50/50 chance of being -ve or +ve
      double p;
      p = mtrandom::grnd();
      if (p >= 0.5)
         factorx *= -1;
      p = mtrandom::grnd();
      if (p >= 0.5)
         factory *= -1;
      p = mtrandom::grnd();
      if (p >= 0.5)
         factorz *= -1;

      new_spin[0] = old_spin[0] + factorx * angle;
      new_spin[1] = old_spin[1] + factory * angle;
      new_spin[2] = old_spin[2] + factorz * angle;

      return;
   }

   int lsf_mc_step(){

      // Enable calling of routine if error checking is activated
      if (err::check == true){
         std::cout << "sim::lsf_mc has been called" << std::endl;
      }

      const int num_atoms = atoms::num_atoms;

      const double kB = 1.3806503e-23;

      // Temporaries
      int atom = 0;
      double Eold = 0.0;
      double Enew = 0.0;
      double DE = 0.0;
      double P = 0.0;

      // Material dependent temperature rescaling
      std::vector<double> rescaled_material_kBTBohr(internal::num_materials);
      std::vector<double> sigma_array(internal::num_materials); // range for tuned gaussian random move
      for (int m = 0; m < internal::num_materials; ++m){
         double alpha = internal::temperature_rescaling_alpha[m];
         double Tc = internal::temperature_rescaling_Tc[m];
         double rescaled_temperature = sim::temperature < Tc ? Tc * pow(sim::temperature / Tc, alpha) : sim::temperature;
         rescaled_material_kBTBohr[m] = 9.27400915e-24 / (rescaled_temperature * 1.3806503e-23);
         sigma_array[m] = rescaled_temperature < 1.0 ? 0.02 : pow(1.0 / rescaled_material_kBTBohr[m], 0.2) * 0.08;
      }

      double statistics_moves = 0.0;
      double statistics_reject = 0.0;

      // Initialise LSF-MC
      if (montecarlo::mc_set == false)
         montecarlo::mcinit();

      for (int i = 0; i < atoms::num_atoms; i++){

         atom = int(atoms::num_atoms * mtrandom::grnd());

         // add one to number of moves counter
         statistics_moves += 1.0;

         // get material id
         const int imaterial = atoms::type_array[atom];

         // Save old spin position
         internal::Sold[0] = atoms::x_spin_array[atom];
         internal::Sold[1] = atoms::y_spin_array[atom];
         internal::Sold[2] = atoms::z_spin_array[atom];

         // Transverse step
         montecarlo::mc_angle_transverse(internal::Sold, internal::Snew, montecarlo::internal::adaptive_sigma, atom);

         atoms::x_spin_array[atom] = internal::Sold[0];
         atoms::y_spin_array[atom] = internal::Sold[1];
         atoms::z_spin_array[atom] = internal::Sold[2];

         // Calculate current energy
         Eold = sim::calculate_spin_energy(atom);

         atoms::x_spin_array[atom] = internal::Snew[0];
         atoms::y_spin_array[atom] = internal::Snew[1];
         atoms::z_spin_array[atom] = internal::Snew[2];

         // Calculate new energy
         Enew = sim::calculate_spin_energy(atom);

         // Calculate difference in Joules/mu_B
         DE = (Enew - Eold) * internal::mu_s_SI[imaterial] * 1.07828231e23; // 1/9.27400915e-24

         double P = exp(-DE * rescaled_material_kBTBohr[imaterial]);

         if (DE < 0){

            // Copy new spin position
            atoms::x_spin_array[atom] = internal::Snew[0];
            atoms::y_spin_array[atom] = internal::Snew[1];
            atoms::z_spin_array[atom] = internal::Snew[2];
         } else if (P >= mtrandom::grnd()){

            // Copy new spin position
            atoms::x_spin_array[atom] = internal::Snew[0];
            atoms::y_spin_array[atom] = internal::Snew[1];
            atoms::z_spin_array[atom] = internal::Snew[2];
         } else{

            atoms::x_spin_array[atom] = internal::Sold[0];
            atoms::y_spin_array[atom] = internal::Sold[1];
            atoms::z_spin_array[atom] = internal::Sold[2];

            // add one to rejection counter
            statistics_reject += 1.0;
         }
      }

      // calculate new adaptive step sigma angle
      if (statistics_moves != 0){
         const double last_rejection_rate = statistics_reject / statistics_moves;
         const double factor = 0.5 / last_rejection_rate;
         montecarlo::internal::adaptive_sigma *= factor;
         // check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
         if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5)
            montecarlo::internal::adaptive_sigma = 60.0; // FOR NORMAL RUNS
      }

      // Save statistics to sim namespace variable
      sim::mc_statistics_moves += statistics_moves;
      sim::mc_statistics_reject += statistics_reject;
   }

} // namespace montecarlo
