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
#include "quantum.hpp"
#include "vmpi.hpp"

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

      // Guard against multiple initialization (can happen in MPI when called from different places)
      static bool already_initialized = false;
      if(already_initialized) return;
      already_initialized = true;

      std::cout << "Initializing Quantum Noise Module..." << std::endl;

      //------------------------------------------------------------------------
      // Common initialization - populate material parameters
      //------------------------------------------------------------------------
      for(int m=0; m < mp::num_materials; m++){

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
         material_S0_array.push_back(S0);
         material_inv_sqrt_S0_array.push_back(inv_sqrt_S0);

      }

      //------------------------------------------------------------------------
      // Disable standard thermal field (index 3) - common to all methods
      //------------------------------------------------------------------------
      sim::hamiltonian_simulation_flags[3] = 0;

      //------------------------------------------------------------------------
      // Check that this program is explicitly supported - common to all methods
      //------------------------------------------------------------------------
      uint64_t total_simulation_time = 0;
      if(!supported_program(total_simulation_time)){
         std::cerr << "Error! Selected program " << program::program << " is not currently supported when using the quantum thermostat" << std::endl;
         err::vexit();
      }

      //------------------------------------------------------------------------
      // Print material parameters - common to all methods
      //------------------------------------------------------------------------
      for(int m=0; m < mp::num_materials; m++){
         std::cout << "  Material " << m << " parameters: A = " << material_A_array[m]
                   << ", Gamma = " << material_gamma_array[m]
                   << ", omega0 = " << material_omega0_array[m] << 
                   ", S0 = " << material_S0_array[m] << std::endl;
      }

      //------------------------------------------------------------------------
      // Dispatch to method-specific initialization
      //------------------------------------------------------------------------
      if(internal::llg_method == internal::llg_ho){
         initialize_HO();
      }
      else if(internal::llg_method == internal::llg_fft){
         initialize_FFT();
      }

      return;

   }

   //----------------------------------------------------------------------------
   // FFT-specific initialization
   //----------------------------------------------------------------------------
   void initialize_FFT(){

      using namespace internal;

      std::cout << "  Using FFT-based quantum noise generation" << std::endl;

      // Get total simulation time
      uint64_t total_simulation_time = 0;
      supported_program(total_simulation_time); // get total simulation time

      // How often calculate random fields? 3 components per atom (x,y,z)
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
      // Calculate coarse grid size based on decimation factor
      int n_coarse = (n_fine > 0) ? ((n_fine - 1) / M_decimation + 1) : 0;

      std::cout << "  FFT simulation parameters:" << std::endl;
      std::cout << "    Temperature: " << T << " K" << std::endl;
      std::cout << "    Time step: " << dt_fine << " s" << std::endl;
      std::cout << "    Total fine steps: " << n_fine << std::endl;
      std::cout << "    Interpolation factor M: " << M_decimation << std::endl;
      std::cout << "    Coarse steps: " << n_coarse << std::endl;
      std::cout << "    Window size: " << window_size << std::endl;
      std::cout << "    Noise type: " << noise_type << std::endl;

      // Determine number of atoms (including boundary atoms for MPI)
      #ifdef MPICF
         const int num_atoms_total = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         const int num_atoms_total = atoms::num_atoms;
      #endif

      // Resize arrays (FFT method uses 9 components like HO: S, q, p)
      q_x_array.resize(num_atoms_total, 0.0);
      q_y_array.resize(num_atoms_total, 0.0);
      q_z_array.resize(num_atoms_total, 0.0);
      p_x_array.resize(num_atoms_total, 0.0);
      p_y_array.resize(num_atoms_total, 0.0);
      p_z_array.resize(num_atoms_total, 0.0);

      k1_storage.resize(num_atoms_total, std::vector<double>(9));
      k2_storage.resize(num_atoms_total, std::vector<double>(9));
      k3_storage.resize(num_atoms_total, std::vector<double>(9));
      k4_storage.resize(num_atoms_total, std::vector<double>(9));
      y_pred_storage.resize(num_atoms_total, std::vector<double>(9));
      y_in_storage.resize(num_atoms_total, std::vector<double>(9));

      // Assign unique indices for random fields (use local atom count for sizing)
      assign_unique_indices(n_coarse, num_atoms_total);

      // Initialize noise structures (FFT plans, buffers, PSDs)
      init_noise_structures(n_coarse, realizations, dt_fine, M_decimation, T);

      // Generate initial random fields
      calculate_noise(realizations,
                      n_fine,
                      dt_fine,
                      M_decimation,
                      T,
                      n_coarse,
                      coarse_noise_field);

      return;

   }

   //----------------------------------------------------------------------------
   // Harmonic Oscillator-specific initialization
   //----------------------------------------------------------------------------
   void initialize_HO(){

      using namespace internal;

      std::cout << "  Using Harmonic Oscillator (direct noise generation)" << std::endl;

      // Determine number of atoms (including boundary atoms for MPI)
      #ifdef MPICF
         const int num_atoms_total = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         const int num_atoms_total = atoms::num_atoms;
      #endif

      // Resize arrays
      q_x_array.resize(num_atoms_total, 0.0);
      q_y_array.resize(num_atoms_total, 0.0);
      q_z_array.resize(num_atoms_total, 0.0);
      p_x_array.resize(num_atoms_total, 0.0);
      p_y_array.resize(num_atoms_total, 0.0);
      p_z_array.resize(num_atoms_total, 0.0);

      k1_storage.resize(num_atoms_total, std::vector<double>(9));
      k2_storage.resize(num_atoms_total, std::vector<double>(9));
      k3_storage.resize(num_atoms_total, std::vector<double>(9));
      k4_storage.resize(num_atoms_total, std::vector<double>(9));
      y_pred_storage.resize(num_atoms_total, std::vector<double>(9));
      y_in_storage.resize(num_atoms_total, std::vector<double>(9));

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
      programs[6] = true;
      simulation_time[6] = et+tt;

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
