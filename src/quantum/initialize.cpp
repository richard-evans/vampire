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

      std::cout << "Initializing Quantum Noise Module" << std::endl;

      // Populate Lorentzian parameter arrays
      for(int i=0; i < atoms::num_atoms; i++){

         // Get material index
         const int imaterial=atoms::type_array[i];

         // Get material parameters
         double alpha = mp::material[imaterial].alpha;
         double gamma = mp::material[imaterial].gamma.get();
         double omega0 = mp::material[imaterial].omega0.get();

         // Calculate amplitude from damping parameter and Gamma and omega0
         double A = alpha * pow(omega0, 4) / gamma;

         // Populate parameter arrays
         material_gamma_array.push_back(gamma);
         material_omega0_array.push_back(omega0);
         material_A_array.push_back(A);
         
      }

   
      // Disable standard thermal field (index 3)
      sim::hamiltonian_simulation_flags[3] = 0;

      // How often calcaulte random fields? 3 components per atom (x,y,z)
      int realizations = atoms::num_atoms * 3;

      // Temperature
      double T = sim::temperature;

      // if window size is zero, set to simulation length
      if(window_size == 0){
         window_size = static_cast<int>(sim::total_time) + 1;
      }

      // Fine time array (used for spin dynamics)
      double dt_fine = mp::dt;
      int n_fine = static_cast<int>(sim::total_time) + 1;

      // Coarse time array (used for noise generation)
      int M_decimation = internal::M_decimation;
      int n_coarse = window_size;
      
      
      std::cout << "Quantum noise module simulation parameters:" << std::endl;
      std::cout << "  Temperature: " << T << " K" << std::endl;
      std::cout << "  Time step: " << dt_fine << " s" << std::endl;
      std::cout << "  Total fine steps: " << n_fine << std::endl;
      std::cout << "  Decimation factor M: " << M_decimation << std::endl;
      std::cout << "  Coarse steps: " << n_coarse << std::endl;
      std::cout << "  Window size: " << window_size << std::endl;
      std::cout << "  Noise type: " << noise_type << std::endl;

      // Check if we have valid parameters (print out parameters for each material)
      for(int i = 0; i < material_A_array.size(); i++){
         std::cout << " Material " << i << " parameters: A = " << material_A_array[i]
                   << ", Gamma = " << material_gamma_array[i]
                   << ", omega0 = " << material_omega0_array[i] << std::endl;

         if(material_A_array[i] <= 0.0 || material_gamma_array[i] <= 0.0 || material_omega0_array[i] <= 0.0){
            std::cerr << "Error: Invalid material parameters for material " << i << std::endl;
            err::vexit();
         }
      }

      // Assign unique indices for random fields
      assign_unique_indices(n_coarse);

      // Generate random fields
      calculate_noise(realizations, 
                      n_fine, 
                      dt_fine, 
                      M_decimation, 
                      T, 
                      n_coarse, 
                      coarse_noise_field);

      return;

   }

} // end of quantum namespace
