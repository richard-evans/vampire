//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "quantum.hpp"
#include "sim.hpp"
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   //----------------------------------------------------------------------------
   // Function to initialize quantum module
   //----------------------------------------------------------------------------
   void initialize(){

      using namespace internal;

      if(!enabled && mp.empty()) return;

      // If mp is empty but enabled is true, create default mp
      if(mp.empty()){
          mp.resize(1);
      }

      // Check if we have valid parameters
      if(!mp[0].A.is_set() || !mp[0].Gamma.is_set() || !mp[0].omega0.is_set()){
          if(enabled){
             std::cout << "Warning: Quantum module enabled but material parameters (A, Gamma, omega0) not fully set." << std::endl;
          }
          return;
      }

      std::cout << "Initializing Quantum Noise Module..." << std::endl;

      // Disable standard thermal field (index 3 usually)
      // sim::hamiltonian_simulation_flags is a vector<int> or array?
      // In llg_quantum.cpp: sim::hamiltonian_simulation_flags[3] = 0;
      sim::hamiltonian_simulation_flags[3] = 0;

      const int num_atoms = atoms::num_atoms;
      // 3 components per atom (x,y,z)
      int realizations = num_atoms * 3;

      noise_index = 0.0;

      // Time step and duration
      // Assuming sim::dt and sim::total_steps are available
      // If not, we might need to include program.hpp or similar
      double dt = ::dt;
      long int n_fine = sim::total_steps;

      // Fallback if total_steps is not set (e.g. equilibration only)
      if(n_fine <= 0) n_fine = static_cast<long int>(sim::equilibration_time) + static_cast<long int>(sim::loop_time);
      if(n_fine <= 0) n_fine = 10000; // Default fallback

      // Temperature
      double T = sim::temperature;

      double omega_cutoff = estimate_cutoff_omega_cdf(T, 0.99999);
      M_decimation = static_cast<int>(std::ceil((M_PI / omega_cutoff) / dt));
      if(M_decimation < 1) M_decimation = 1;

      int n_coarse = (n_fine > 0) ? ((n_fine - 1) / M_decimation + 1) : 0;

      std::cout << "  Quantum Noise Configuration:" << std::endl;
      std::cout << "  Temperature: " << T << " K" << std::endl;
      std::cout << "  Time step: " << dt << " s" << std::endl;
      std::cout << "  Total fine steps: " << n_fine << std::endl;
      std::cout << "  Decimation factor M: " << M_decimation << std::endl;
      std::cout << "  Coarse steps: " << n_coarse << std::endl;
      std::cout << "  Window size: " << window_size << std::endl;

      // Assign unique indices for random fields
      assign_unique_indices(n_coarse);

      // Generate random fields
      calculate_random_fields_windowed(realizations, n_fine, dt, M_decimation, T, n_coarse, coarse_noise_field);

      return;

   }

} // end of quantum namespace
