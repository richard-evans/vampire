//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "vutil.hpp"

// Include internal header
#include "internal.hpp"

namespace montecarlo{

//------------------------------------------------------------------------------
// Function to perform simple Monte Carlo solver as a preconditioner to
// obtain a thermodynamically sensible initial state for the spins at the
// equilibration temperature. This is typically > 10000 times more efficient
// than direct LLG integration.
//
// Note: this algorithm does not maintain ergodicity in parallel and is only
// definitively safe to use as a preconditioner rather than for full dynamics.
//------------------------------------------------------------------------------
void monte_carlo_preconditioning(){

   //if no steps specified, then do nothing
   if(sim::num_monte_carlo_preconditioning_steps == 0) return;

   // instantiate timer
   vutil::vtimer_t timer;

   // start timer
   timer.start();

   // set preconditioning temperature to equilibration temperature
   sim::temperature = sim::Teq;

   // print informative messages to screen and log
   std::cout << "Preconditioning spin configuration at T = " << sim::temperature << " K" << std::flush;
   zlog << zTs() << "Preconditioning spin configuration at T = " << sim::temperature << " K ..." << std::endl;

   // calculate number of steps to calculate
   // Note: In parallel this includes boundary spin leading to
   // spurious edge dynamics
	const int num_moves = atoms::num_atoms;

   // Material dependent temperature rescaling unrolled for speed
   std::vector<double> rescaled_material_kBTBohr(mp::num_materials);
   std::vector<double> moment_array(mp::num_materials); // mu_s/mu_B
   std::vector<double> sigma_array(mp::num_materials); // range for tuned gaussian random move
   for(int m=0; m < mp::num_materials; ++m){
      double alpha = mp::material[m].temperature_rescaling_alpha;
      double Tc = mp::material[m].temperature_rescaling_Tc;
      double rescaled_temperature = sim::temperature < Tc ? Tc*pow(sim::temperature/Tc,alpha) : sim::temperature;
      rescaled_material_kBTBohr[m] = 9.27400915e-24/(rescaled_temperature*1.3806503e-23);
      sigma_array[m] = rescaled_temperature < 1.0 ? 0.02 : pow(1.0/rescaled_material_kBTBohr[m],0.2)*0.08;
      moment_array[m] = mp::material[m].mu_s_SI/9.27400915e-24;
   }

	for(int s = 0; s < sim::num_monte_carlo_preconditioning_steps; s++){

      //std::cout << s << std::endl;
      if( (s % (sim::num_monte_carlo_preconditioning_steps/10)) == 0) std::cout << "." << std::flush;

      // loop over natoms to form a single Monte Carlo step
      for(int i = 0; i < num_moves; i++){

   		// pick atom
   		const int atom = int(num_moves * mtrandom::grnd());

         //std::cout << atom << std::endl;

   		// get material id
   		const int imaterial = atoms::type_array[atom];

         //std::cout << "here-2 " << imaterial << std::endl;

         // Calculate range for move
         internal::delta_angle = sigma_array[imaterial];

         //std::cout << "here-1 " << sim::mc_delta_angle << std::endl;

   		// Save old spin position
   		internal::Sold[0] = atoms::x_spin_array[atom];
   		internal::Sold[1] = atoms::y_spin_array[atom];
   		internal::Sold[2] = atoms::z_spin_array[atom];

         //std::cout << "here0" << std::endl;

         // Make Monte Carlo move
         montecarlo::internal::mc_move(internal::Sold, internal::Snew);

         //std::cout << "here" << std::endl;

   		// Calculate current energy
   		double old_energy = sim::calculate_spin_energy(atom);

   		// Copy new spin position
   		atoms::x_spin_array[atom] = internal::Snew[0];
   		atoms::y_spin_array[atom] = internal::Snew[1];
   		atoms::z_spin_array[atom] = internal::Snew[2];

         //std::cout << "here2" << std::endl;


   		// Calculate new energy
   		const double new_energy = sim::calculate_spin_energy(atom);

   		// Calculate difference in Joules/mu_B
   		const double DE = (new_energy - old_energy) * moment_array[imaterial];

   		// Check for lower energy state and accept unconditionally
   		if(DE<0) continue;
   		// Otherwise evaluate probability for move
   		else{
   			if(exp(-DE*rescaled_material_kBTBohr[imaterial]) >= mtrandom::grnd()) continue;
   			// If rejected reset spin coordinates and continue
   			else{
   				atoms::x_spin_array[atom] = internal::Sold[0];
   				atoms::y_spin_array[atom] = internal::Sold[1];
   				atoms::z_spin_array[atom] = internal::Sold[2];
   				continue;
   			}
   		}

   	} // end of each MC move

      // update halo data
      vmpi::mpi_init_halo_swap();
      vmpi::mpi_complete_halo_swap();

   } // end of preconditioning steps loop

   // start timer
   timer.stop();

   // print informative messages to screen and log
   std::cout << "Done!" << std::endl;
   std::cout << "Preconditioning time for " << sim::num_monte_carlo_preconditioning_steps << " steps: " << timer.elapsed_time() << " s" << std::endl;
   zlog << zTs() << "Preconditioning completed in " << timer.elapsed_time() << " s" << std::endl;

   return;

}

} // end of namespace montecarlo
