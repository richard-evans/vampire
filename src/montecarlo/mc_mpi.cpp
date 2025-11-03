//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans and Adam Laverack 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifdef MPICF
// Standard Libraries
#include <cmath>
#include <iostream>
#include <vector>

// Vampire Header files
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Internal header
#include "internal.hpp"

namespace montecarlo{

//------------------------------------------------------------------------------
// Initialise octant arrays to store which atom is in which octant for the
// checkerboard MC algorithm
//------------------------------------------------------------------------------
void mc_parallel_init(std::vector<double> &x, // atomic coordinates
                      std::vector<double> &y,
                      std::vector<double> &z,
                      double min_dim[3], // minimum dimensions on local processor
                      double max_dim[3]){ // maximum dimensions on local processor

   // Convenient shorthands
   int catoms = vmpi::num_core_atoms;
   int batoms = vmpi::num_bdry_atoms;

   double widthx = max_dim[0] - min_dim[0];
   double widthy = max_dim[1] - min_dim[1];
   double widthz = max_dim[2] - min_dim[2];

   int octant_num = 0; //Count which octant loop is in

   // Determines which core atoms are in which octant and pushes the index of those
   // atoms into the appropriate octant arrays.
   for(int zoct=0; zoct<2; zoct++){
      for(int yoct=0; yoct<2; yoct++){
         for(int xoct=0; xoct<2; xoct++){
            // Loop through all core atoms
            for (int i=0; i<catoms; i++){
               // Check if current atom is in desired octant
               if (   x[i] >= min_dim[0] + widthx*xoct*0.5 && x[i] < min_dim[0] + widthx*0.5 + widthx*xoct*0.5
                   && y[i] >= min_dim[1] + widthy*yoct*0.5 && y[i] < min_dim[1] + widthy*0.5 + widthy*yoct*0.5
                   && z[i] >= min_dim[2] + widthz*zoct*0.5 && z[i] < min_dim[2] + widthz*0.5 + widthz*zoct*0.5)
               {
                  internal::c_octants[octant_num].push_back(i);
               }
            }
            octant_num++;
         }
      }
   }

   octant_num = 0;
   //Sort boundary atoms into appropriate octant arrays.
   for(int zoct=0; zoct<2; zoct++){
      for(int yoct=0; yoct<2; yoct++){
         for(int xoct=0; xoct<2; xoct++){
            // Loop through all boundary atoms
            for (int i=catoms; i<catoms+batoms; i++){
               // Check if current atom is in desired octant
               if (   x[i] >= min_dim[0] + widthx*xoct*0.5 && x[i] < min_dim[0] + widthx*0.5 + widthx*xoct*0.5
                   && y[i] >= min_dim[1] + widthy*yoct*0.5 && y[i] < min_dim[1] + widthy*0.5 + widthy*yoct*0.5
                   && z[i] >= min_dim[2] + widthz*zoct*0.5 && z[i] < min_dim[2] + widthz*0.5 + widthz*zoct*0.5)
               {
                  internal::b_octants[octant_num].push_back(i);
               }
            }
            octant_num++;
         }
      }
   }

   //--------------------------------------------------------------------
   // check that all atoms have been allocated an octant
   //--------------------------------------------------------------------
   // core atoms
   int num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) num_atoms_in_octants += internal::c_octants[i].size();
   if(num_atoms_in_octants != catoms){
      std::cerr << "Programmer error: missing atoms in core octants in parallel monte carlo initialisation" << std::endl;
      err::vexit();
   }
   // boundary atoms
   num_atoms_in_octants = 0;
   for(int i=0; i< 8; i++) num_atoms_in_octants += internal::b_octants[i].size();
   if(num_atoms_in_octants != batoms){
      std::cerr << "Programmer error: missing atoms in boundary octants in parallel monte carlo initialisation" << std::endl;
      err::vexit();
   }

   // set flag to indicate that parallel MC has been initialised
   mc_parallel_initialized = true;

}

//------------------------------------------------------------------------------
// Integrates a Monte Carlo step in parallel
//------------------------------------------------------------------------------
void mc_step_parallel(std::vector<double> &x_spin_array,
                      std::vector<double> &y_spin_array,
                      std::vector<double> &z_spin_array,
                      std::vector<int> &type_array){

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

	// loop over all octants
   for(int octant = 0; octant < 8; octant++) {

      int nmoves = internal::c_octants[octant].size();
      //Initialise non-blocking send
      vmpi::mpi_init_halo_swap();

      //---------------------------------
      // Integrate core region
      //---------------------------------
      for(int i=0; i<nmoves; i++){

         // add one to number of moves counter
         statistics_moves+=1.0;

      	// pick random atom in octant
      	atom = internal::c_octants[octant][int(nmoves*mtrandom::grnd())];

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

      // Finish non-blocking data send/receive
      vmpi::mpi_complete_halo_swap();

      //
      //Begin integrating boundary region
      //
      nmoves = internal::b_octants[octant].size();
      for(int i=0; i<nmoves; i++){

         // add one to number of moves counter
         statistics_moves+=1.0;

   		// pick atom
   		atom = internal::b_octants[octant][int(nmoves*mtrandom::grnd())];

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

      // Swap timers compute -> wait
      vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

      // Wait for other processors
      vmpi::barrier();

      // Swap timers wait -> compute
      vmpi::TotalWaitTime += vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

   } // end of octant loop

   //Collect statistics from all processors
   double global_statistics_moves = 0.0;
   double global_statistics_reject = 0.0;
   MPI_Allreduce(&statistics_moves, &global_statistics_moves, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&statistics_reject, &global_statistics_reject, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   // calculate new adaptive step sigma angle (on per-processor basis using local, not global stats)
   if(montecarlo::algorithm == montecarlo::adaptive){
      const double last_rejection_rate = statistics_reject / statistics_moves;
      const double factor = 0.5 / last_rejection_rate;
      montecarlo::internal::adaptive_sigma *= factor;
      // check for excessive range (too small angle takes too long to grow, too large does not improve performance) and truncate
      if (montecarlo::internal::adaptive_sigma > 60.0 || montecarlo::internal::adaptive_sigma < 1e-5) montecarlo::internal::adaptive_sigma = 60.0;
   }

   // Save statistics to sim namespace variable
   sim::mc_statistics_moves += global_statistics_moves;
   sim::mc_statistics_reject += global_statistics_reject;

   return;

}

} // End of namespace montecarlo
#endif
