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
#include "program.hpp"

// quantum module headers
#include "internal.hpp"

namespace quantum{

   bool supported_program(uint64_t& total_simulation_time);
   void initialize_FFT(); // Forward declaration for FFT-specific initialization
   void initialize_HO();  // Forward declaration for HO-specific initialization

   //----------------------------------------------------------------------------
   // Common initialization wrapper for quantum module
   //----------------------------------------------------------------------------
   void initialize(){

      using namespace internal;

      // Do nothing if module not enabled
      if(!enabled) return;

      std::cout << "Initializing Quantum Noise Module" << std::endl;

      // Populate Lorentzian parameter arrays
      for(int m=0; m < mp::num_materials; m++){

         // Get material index
         //const int imaterial=atoms::type_array[i];

         // Get material parameters
         double alpha = mp::material[m].alpha; // get damping constant from main materials class
         double gamma = internal::mp[m].gamma.get();
         double omega0 = internal::mp[m].omega0.get();

         // Get material parameters
         double alpha = mp::material[m].alpha; // get damping constant from main materials class
         double gamma = internal::mp[m].gamma.get();
         double omega0 = internal::mp[m].omega0.get();
         double S0 = mp::material[m].mu_s_SI / 9.274009994e-24; // spin magnitude from mu_s and mu_B
         double inv_sqrt_S0 = (S0 > 0) ? 1.0 / std::sqrt(S0) : 1.0;
         // Calculate amplitude from damping parameter and Gamma and omega0
         double A = alpha * pow(omega0, 4) / gamma;

         // Populate parameter arrays
         material_gamma_array.push_back(gamma);
         material_omega0_array.push_back(omega0);
         material_A_array.push_back(A);

      }


      // Disable standard thermal field (index 3)
      sim::hamiltonian_simulation_flags[3] = 0;


      // Setting up the window
      uint64_t total_simulation_time = 0;

      // Check that this program is explicitly supported
      if(!supported_program(total_simulation_time)){
         // error
         std::cerr << "Error! Selected program " << program::program << "is not currently supported when using the quantum thermostat" << std::endl;
         err::vexit();
      }

      // How often calcaulte random fields? 3 components per atom (x,y,z)
      int realizations = atoms::num_atoms * 3;

      // Temperature
      double T = sim::temperature;

      // if window size is zero, set to simulation length
      if(window_size == 0){
         window_size = total_simulation_time + 1;
      }

      // Fine time array (used for spin dynamics)
      double dt_fine = mp::dt;
      uint64_t n_fine = sim::total_time + 1;

      // Coarse time array (used for noise generation)
      int M_decimation = internal::M_decimation;
      int n_coarse = window_size;

      std::cout << "Quantum noise module simulation parameters:" << std::endl;
      std::cout << "  Temperature: " << T << " K" << std::endl;
      std::cout << "  Time step: " << dt_fine << " s" << std::endl;
      std::cout << "  Total fine steps: " << n_fine << std::endl;
      std::cout << "  Interpolation factor M: " << M_decimation << std::endl;
      std::cout << "  Coarse steps: " << n_coarse << std::endl;
      std::cout << "  Window size: " << window_size << std::endl;
      std::cout << "  Noise type: " << noise_type << std::endl;

      // Check if we have valid parameters (print out parameters for each material)
      for(int m=0; m < mp::num_materials; m++){
         std::cout << " Material " << m << " parameters: A = " << material_A_array[m]
                   << ", Gamma = " << material_gamma_array[m]
                   << ", omega0 = " << material_omega0_array[m] << std::endl;

         //if(material_A_array[m] <= 0.0 || material_gamma_array[m] <= 0.0 || material_omega0_array[m] <= 0.0){
         //   std::cerr << "Error: Invalid material parameters for material " << i << std::endl;
         //   err::vexit();
         //}
      }

      // Resize arrays
      x_w_array.resize(atoms::num_atoms, 0.0);
      y_w_array.resize(atoms::num_atoms, 0.0);
      z_w_array.resize(atoms::num_atoms, 0.0);
      x_v_array.resize(atoms::num_atoms, 0.0);
      y_v_array.resize(atoms::num_atoms, 0.0);
      z_v_array.resize(atoms::num_atoms, 0.0);

      k1_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k2_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k3_storage.resize(atoms::num_atoms, std::vector<double>(9));
      k4_storage.resize(atoms::num_atoms, std::vector<double>(9));
      y_pred_storage.resize(atoms::num_atoms, std::vector<double>(9));
      y_in_storage.resize(atoms::num_atoms, std::vector<double>(9));

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

      // Set module as initialised
      internal::initialised = true;

      return;

   }

   bool supported_program(uint64_t& total_simulation_time){

      // set a list of programs assuming a maximum large number
      std::vector<bool> programs(1000, false);
      std::vector<uint64_t> simulation_time(1000, 0);

      if(program::program >= 1000){
         std::cerr << "Programmer error! Program ID in quantum module initialisation exceeds the maximum number of 1000. Decrease program ID or increase number of supported programs" << std::endl;
         err::vexit();
      }

      const uint64_t et = sim::equilibration_time;
      const uint64_t tt = sim::total_time;
      const uint64_t lt = sim::loop_time;


      // Benchmark (total time)
      programs[0] = true;
      simulation_time[0] = tt;
      //programs[program::benchmark] = true;

      // Time series (eq+total time)
      programs[1] = true;
      simulation_time[1] = et+tt;

      // Hysteresis (eq+loop time)
      programs[2] = true;
      simulation_time[2] = et+lt;

      // Static hysteresis (eq+loop time)
      programs[3] = true;
      simulation_time[3] = et+lt;

      // Curie temperature (eq+loop time)
      programs[4] = true;
      simulation_time[4] = et+lt;

      // Field cool (not supported - dynamic temperature)
      programs[5] = false;
      //simulation_time[5] = et+lt;

      // Temperature pulse (not supported, dynamic temperature)
      programs[6] = false;

      // HAMR (not supported, dynamic temperature)
      programs[7] = false;

      // LaGrange multiplier (eq+total time)
      programs[11] = true;
      simulation_time[11] = et+tt;

      // Partial hysteresis (eq+loop time)
      programs[12] = true;
      simulation_time[12] = et+lt;

      // Localised temperature pulse (not supported, dynamic temperature)
      programs[13] = false;

      // Effective damping (eq+total time)
      programs[14] = true;
      simulation_time[14] = et+tt;

      // FMR (eq+total time)
      programs[15] = true;
      simulation_time[15] = et+tt;

      // Local Field cool (not supported - dynamic temperature)
      programs[16] = false;
      //simulation_time[16] = et+lt;

      // Electrical pulse (eq+total time)
      programs[17] = true;
      simulation_time[17] = et+tt;

      // Field pulse (eq+total time)
      programs[18] = true;
      simulation_time[18] = et+tt;

      // Domain walls (eq+total time)
      programs[52] = true;
      simulation_time[52] = et+tt;

      // Field sweep (eq+loop time)
      programs[70] = true;
      simulation_time[70] = et+lt;

      // Spin waves (eq+total time)
      programs[74] = true;
      simulation_time[74] = et+tt;

      total_simulation_time = simulation_time[program::program];

      // return a bool that is true if the current selected program is supported by the quantum thermostat
      return programs[program::program];

   }

} // end of quantum namespace
