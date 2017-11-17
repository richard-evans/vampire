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
#include <cstdlib>
#include <iostream>
#include <vector>

// Vampire Header files
#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// Internal header
#include "internal.hpp"

namespace montecarlo{

namespace internal{

//------------------------------------------------------------------------------
// Initialise octant arrays to store which atom is in which octant for the
// checkerboard MC algorithm
//------------------------------------------------------------------------------
void mc_parallel_init(){
   // Convenient shorthands
   int &catoms = vmpi::num_core_atoms;
   int &batoms = vmpi::num_bdry_atoms;
   std::vector<double> &x = atoms::x_coord_array;
   std::vector<double> &y = atoms::y_coord_array;
   std::vector<double> &z = atoms::z_coord_array;

   double *min_dim = vmpi::min_dimensions;
   double *max_dim = vmpi::max_dimensions;
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
               if (   x[i] > min_dim[0] + widthx/2.0*xoct && x[i] < widthx/2.0 + widthx/2.0*xoct
                   && y[i] > min_dim[1] + widthy/2.0*yoct && y[i] < widthy/2.0 + widthy/2.0*yoct
                   && z[i] > min_dim[2] + widthz/2.0*zoct && z[i] < widthz/2.0 + widthz/2.0*zoct)
               {
                  c_octants[octant_num].push_back(i);
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
               if (   x[i] > min_dim[0] + widthx/2.0*xoct && x[i] < widthx/2.0 + widthx/2.0*xoct
                   && y[i] > min_dim[1] + widthy/2.0*yoct && y[i] < widthy/2.0 + widthy/2.0*yoct
                   && z[i] > min_dim[2] + widthz/2.0*zoct && z[i] < widthz/2.0 + widthz/2.0*zoct)
               {
                  b_octants[octant_num].push_back(i);
               }
            }
            octant_num++;
         }
      }
   }
   mc_parallel_initialized = true;
}

} //end of namespace internal

   //------------------------------------------------------------------------------
   // Integrates a Monte Carlo step in parallel
   //------------------------------------------------------------------------------
int mc_step_parallel(){

	// Check for calling of function
	if(err::check==true) std::cout << "montecarlo::mc_step has been called" << std::endl;

   if(internal::mc_parallel_initialized == false) {
      internal::mc_parallel_init();
   }

	// Declare arrays for spin states
	std::valarray<double> Sold(3);
	std::valarray<double> Snew(3);

	// Temporaries
	int atom=0;
	double Eold=0.0;
	double Enew=0.0;
	double DE=0.0;

   // Material dependent temperature rescaling
   std::vector<double> rescaled_material_kBTBohr(mp::num_materials);
   std::vector<double> sigma_array(mp::num_materials); // range for tuned gaussian random move
   for(int m=0; m<mp::num_materials; ++m){
      double alpha = mp::material[m].temperature_rescaling_alpha;
      double Tc = mp::material[m].temperature_rescaling_Tc;
      double rescaled_temperature = sim::temperature < Tc ? Tc*pow(sim::temperature/Tc,alpha) : sim::temperature;
      rescaled_material_kBTBohr[m] = 9.27400915e-24/(rescaled_temperature*1.3806503e-23);
      sigma_array[m] = rescaled_temperature < 1.0 ? 0.02 : pow(1.0/rescaled_material_kBTBohr[m],0.2)*0.08;
   }

   double statistics_moves = 0.0;
   double statistics_reject = 0.0;

	// loop over all octants
   for(int octant = 0; octant < 8; octant++) {
      int nmoves = internal::c_octants[octant].size();
      vmpi::mpi_init_halo_swap();
      //loop over core atoms in current octant to begin single monte carlo step
   	for(int i=0; i<nmoves; i++){

         // add one to number of moves counter
         statistics_moves+=1.0;

   		// pick atom
   		atom = internal::c_octants[octant][int(nmoves*mtrandom::grnd())];

   		// get material id
   		const int imaterial=atoms::type_array[atom];

         // Calculate range for move
         sim::mc_delta_angle=sigma_array[imaterial];

   		// Save old spin position
   		Sold[0] = atoms::x_spin_array[atom];
   		Sold[1] = atoms::y_spin_array[atom];
   		Sold[2] = atoms::z_spin_array[atom];

         // Make Monte Carlo move
         internal::mc_move(Sold, Snew);

   		// Calculate current energy
   		Eold = sim::calculate_spin_energy(atom);

   		// Copy new spin position
   		atoms::x_spin_array[atom] = Snew[0];
   		atoms::y_spin_array[atom] = Snew[1];
   		atoms::z_spin_array[atom] = Snew[2];

   		// Calculate new energy
   		Enew = sim::calculate_spin_energy(atom);

   		// Calculate difference in Joules/mu_B
   		DE = (Enew-Eold)*mp::material[imaterial].mu_s_SI*1.07828231e23; //1/9.27400915e-24

   		// Check for lower energy state and accept unconditionally
   		if(DE<0) continue;
   		// Otherwise evaluate probability for move
   		else{
   			if(exp(-DE*rescaled_material_kBTBohr[imaterial]) >= mtrandom::grnd()) continue;
   			// If rejected reset spin coordinates and continue
   			else{
   				atoms::x_spin_array[atom] = Sold[0];
   				atoms::y_spin_array[atom] = Sold[1];
   				atoms::z_spin_array[atom] = Sold[2];
               // add one to rejection counter
               statistics_reject += 1.0;
   				continue;
   			}
   		}
   	}

      vmpi::mpi_complete_halo_swap();

      nmoves = internal::b_octants[octant].size();

      //Loop over all atoms in boundary region to complete monte carlo step
      for(int i=0; i<nmoves; i++){

         // add one to number of moves counter
         statistics_moves+=1.0;

   		// pick atom
   		atom = internal::b_octants[octant][int(nmoves*mtrandom::grnd())];

   		// get material id
   		const int imaterial=atoms::type_array[atom];

         // Calculate range for move
         sim::mc_delta_angle=sigma_array[imaterial];

   		// Save old spin position
   		Sold[0] = atoms::x_spin_array[atom];
   		Sold[1] = atoms::y_spin_array[atom];
   		Sold[2] = atoms::z_spin_array[atom];

         // Make Monte Carlo move
         internal::mc_move(Sold, Snew);

   		// Calculate current energy
   		Eold = sim::calculate_spin_energy(atom);

   		// Copy new spin position
   		atoms::x_spin_array[atom] = Snew[0];
   		atoms::y_spin_array[atom] = Snew[1];
   		atoms::z_spin_array[atom] = Snew[2];

   		// Calculate new energy
   		Enew = sim::calculate_spin_energy(atom);

   		// Calculate difference in Joules/mu_B
   		DE = (Enew-Eold)*mp::material[imaterial].mu_s_SI*1.07828231e23; //1/9.27400915e-24

   		// Check for lower energy state and accept unconditionally
   		if(DE<0) continue;
   		// Otherwise evaluate probability for move
   		else{
   			if(exp(-DE*rescaled_material_kBTBohr[imaterial]) >= mtrandom::grnd()) continue;
   			// If rejected reset spin coordinates and continue
   			else{
   				atoms::x_spin_array[atom] = Sold[0];
   				atoms::y_spin_array[atom] = Sold[1];
   				atoms::z_spin_array[atom] = Sold[2];
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
      vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);
   }

   // Save statistics to sim namespace variable
   sim::mc_statistics_moves += statistics_moves;
   sim::mc_statistics_reject += statistics_reject;



	return EXIT_SUCCESS;
}

} // End of namespace montecarlo
#endif