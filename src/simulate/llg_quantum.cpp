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

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <complex>
#include <iomanip>
#include <numeric>

// Library for FFT
#ifdef FFT
#include <fftw3.h>
#endif

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"

#include "internal.hpp"



namespace sim{

   /// @brief LLG Initialisation function
   ///
   /// @details Resizes arrays used for Heun integration
   ///
   /// @section License
   /// Use of this code, either in source or compiled form, is subject to license from the authors.
   /// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
   ///
   /// @section Information
   /// @author  Richard Evans, richard.evans@york.ac.uk
   /// @version 1.0
   /// @date    07/02/2011
   ///
   /// @return EXIT_SUCCESS
   ///
   /// @internal
   ///	Created:		07/02/2011
   ///	Revision:	  ---
   ///=====================================================================================
   ///

   //---------------------------------------------------------------------------
   // Class for FFTW plan management
   //---------------------------------------------------------------------------
   void calculate_random_fields(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse);
   void assign_unique_indices(int n_coarse);
   void precompute_sqrt_PSD(int n, double dt, double T);
   double PSD(const double& omega, const double& T);
   double estimate_cutoff_omega_cdf(double T, double target_frac);
   double get_noise(const std::vector<double>& coarse_noise, double fine_step_idx, int M, size_t atom_idx);



   int LLGQinit(){
      // check calling of routine if error checking is activated
      if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}

      using namespace LLGQ_arrays;

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


      // Set number of realizations (full field, for final release allow even smaller number of realizations)
      const int num_atoms = atoms::num_atoms;
      int realizations = num_atoms * 3 + 4;
      LLGQ_arrays::noise_index = 0;

      // Disable external thermal field calculations
      sim::hamiltonian_simulation_flags[3] = 0;

      // --- Interpolation Setup ---
      const double dt_fine = mp::dt;
      const int n_fine = static_cast<int>(sim::equilibration_time) + 1;

      // Calculate cutoff frequency and decimation factor
      double omega_cutoff = estimate_cutoff_omega_cdf(sim::temperature, 0.99999);
      M_decimation = static_cast<int>(std::ceil((M_PI / omega_cutoff) / dt_fine));

      const int n_coarse = (n_fine > 0) ? ((n_fine - 1) / M_decimation + 1) : 0;
      double mem_red = 100.0 * (1.0 - static_cast<double>(n_coarse) / n_fine);
      std::cout << "Quantum noise interpolation enabled." << std::endl;
      std::cout << "Decimation factor M=" << M_decimation << ", estimated memory reduction=" << std::fixed << std::setprecision(1) << mem_red << "%" << std::endl;

      // Assign unique indices for random fields
      assign_unique_indices(n_coarse);

      // Generate random fields directly on coarse grid and store for interpolation
      calculate_random_fields(realizations, n_fine, dt_fine, M_decimation, sim::temperature, n_coarse);

      LLG_set=true;

      return EXIT_SUCCESS;
   }

   namespace internal{

    

      /// @brief LLG Heun Integrator Corrector
      ///
      /// @details Integrates the system using the LLG and Heun solver
      // ... (Doxygen comments remain the same) ...
      //
      void spinDynamics(const double* y,const double* H,double* dydt);

      void llg_quantum_step(){

         // check calling of routine if error checking is activated
         if(err::check==true){std::cout << "sim::mLLG has been called" << std::endl;}

         using namespace LLGQ_arrays;

         // Check for initialisation of LLG integration arrays
         if(LLG_set==false) sim::LLGQinit();

         // Local variables
         const int num_atoms = atoms::num_atoms;
         std::vector<double> H(3);

         // Store initial state
         for (int atom = 0; atom < num_atoms; atom++) {
            y_in_storage[atom][0] = atoms::x_spin_array[atom];
            y_in_storage[atom][1] = atoms::y_spin_array[atom];
            y_in_storage[atom][2] = atoms::z_spin_array[atom];
            y_in_storage[atom][3] = LLGQ_arrays::x_v_array[atom];
            y_in_storage[atom][4] = LLGQ_arrays::y_v_array[atom];
            y_in_storage[atom][5] = LLGQ_arrays::z_v_array[atom];
            y_in_storage[atom][6] = LLGQ_arrays::x_w_array[atom];
            y_in_storage[atom][7] = LLGQ_arrays::y_w_array[atom];
            y_in_storage[atom][8] = LLGQ_arrays::z_w_array[atom];
         }


         // Calculate fields
         calculate_spin_fields(0, num_atoms);
         calculate_external_fields(0, num_atoms);

         const double dt = mp::dt;
         const double half_dt = 0.5 * mp::dt;
         const double dt_over_6 = mp::dt / 6.0;
         const int M = M_decimation;

         // K1 Step
         for (int atom = 0; atom < num_atoms; ++atom) {
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index, M, atom_idx_x[atom]);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index, M, atom_idx_y[atom]);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index, M, atom_idx_z[atom]);

            // Calculate K1
            spinDynamics(y_in_storage[atom].data(), H.data(), k1_storage[atom].data());

            // Calculate y_pred for k2
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k1_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;

            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];
         }

         // Update spin fields for k2
         calculate_spin_fields(0, num_atoms);

         // K2 Step
         for (int atom = 0; atom < num_atoms; ++atom) {
            // Update fields for k2 with interpolated noise at time t + dt/2
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_x[atom]);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_y[atom]);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_z[atom]);

            // Calculate K2
            spinDynamics(y_pred_storage[atom].data(), H.data(), k2_storage[atom].data());

            // Calculate y_pred for k3
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k2_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;

            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];

         }

         // Update spin fields for k3
         calculate_spin_fields(0, num_atoms);

         // K3 Step
         for (int atom = 0; atom < num_atoms; ++atom) {
            // Update fields for k3 with interpolated noise at time t + dt/2
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_x[atom]);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_y[atom]);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index+0.5, M, atom_idx_z[atom]);

            // Calculate K3
            spinDynamics(y_pred_storage[atom].data(), H.data(), k3_storage[atom].data());

            // Calculate y_pred for k4
            for (size_t i = 0; i < 9; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + dt * k3_storage[atom][i];
            }

            // Normalize spin length
            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;

            // Update spin for field calculation
            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];
         }

         // Update fields for k4
         calculate_spin_fields(0, num_atoms);

         //K4 step
         for (int atom = 0; atom < num_atoms; ++atom) {
            // Update fields for k4 with interpolated noise at t+dt
            H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index + 1.0, M, atom_idx_x[atom]);
            H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index + 1.0, M, atom_idx_y[atom]);
            H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + get_noise(coarse_noise_field, noise_index + 1.0, M, atom_idx_z[atom]);

            // Calculate K4
            spinDynamics(y_pred_storage[atom].data(), H.data(), k4_storage[atom].data());

            // Final update and normalization for spin components
            for (size_t i = 0; i < 3; ++i) {
               y_pred_storage[atom][i] = y_in_storage[atom][i] + dt_over_6 * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
            }

            double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
            const double inv_mag = 1.0 / S_magnitude;
            y_pred_storage[atom][0] *= inv_mag;
            y_pred_storage[atom][1] *= inv_mag;
            y_pred_storage[atom][2] *= inv_mag;

            atoms::x_spin_array[atom] = y_pred_storage[atom][0];
            atoms::y_spin_array[atom] = y_pred_storage[atom][1];
            atoms::z_spin_array[atom] = y_pred_storage[atom][2];

            // Update auxiliary variables
            LLGQ_arrays::x_v_array[atom] = y_in_storage[atom][3] + (dt_over_6) * (k1_storage[atom][3] + 2.0 * k2_storage[atom][3] + 2.0 * k3_storage[atom][3] + k4_storage[atom][3]);
            LLGQ_arrays::y_v_array[atom] = y_in_storage[atom][4] + (dt_over_6) * (k1_storage[atom][4] + 2.0 * k2_storage[atom][4] + 2.0 * k3_storage[atom][4] + k4_storage[atom][4]);
            LLGQ_arrays::z_v_array[atom] = y_in_storage[atom][5] + (dt_over_6) * (k1_storage[atom][5] + 2.0 * k2_storage[atom][5] + 2.0 * k3_storage[atom][5] + k4_storage[atom][5]);
            LLGQ_arrays::x_w_array[atom] = y_in_storage[atom][6] + (dt_over_6) * (k1_storage[atom][6] + 2.0 * k2_storage[atom][6] + 2.0 * k3_storage[atom][6] + k4_storage[atom][6]);
            LLGQ_arrays::y_w_array[atom] = y_in_storage[atom][7] + (dt_over_6) * (k1_storage[atom][7] + 2.0 * k2_storage[atom][7] + 2.0 * k3_storage[atom][7] + k4_storage[atom][7]);
            LLGQ_arrays::z_w_array[atom] = y_in_storage[atom][8] + (dt_over_6) * (k1_storage[atom][8] + 2.0 * k2_storage[atom][8] + 2.0 * k3_storage[atom][8] + k4_storage[atom][8]);

         }

         // Increment noise index
         LLGQ_arrays::noise_index += 1;

         return;
      }

      void spinDynamics(const double* y, const double* H, double* dydt) {
         const double A = sim::internal::mp[0].A.get();
         const double Gamma = sim::internal::mp[0].Gamma.get();
         const double omega0 = sim::internal::mp[0].omega0.get();

         // dS/dt = S × (H + v)
         dydt[0] =  (y[1]*(H[2]+y[5]) - y[2]*(H[1]+y[4]));
         dydt[1] =  (y[2]*(H[0]+y[3]) - y[0]*(H[2]+y[5]));
         dydt[2] =  (y[0]*(H[1]+y[4]) - y[1]*(H[0]+y[3]));


         // dv/dt = w
         dydt[3] = y[6];
         dydt[4] = y[7];
         dydt[5] = y[8];

         // dw/dt = -ω₀²v - Γw + AS
         dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0];
         dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1];
         dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2];
      }

   } // end of internal namespace

   void assign_unique_indices(int n_coarse) {
      const int num_atoms = atoms::num_atoms;

      std::cout << "Assigning indices for " << num_atoms << " atoms with " << n_coarse << " coarse steps." << std::endl;

      LLGQ_arrays::atom_idx_x.resize(num_atoms);
      LLGQ_arrays::atom_idx_y.resize(num_atoms);
      LLGQ_arrays::atom_idx_z.resize(num_atoms);

      for (int atom = 0; atom < num_atoms; atom++) {
         // Use 64-bit arithmetic to prevent overflow
         const size_t atom_ll = static_cast<size_t>(atom);
         const size_t n_coarse_ll = static_cast<size_t>(n_coarse);
         LLGQ_arrays::atom_idx_x[atom] = 3 * atom_ll * n_coarse_ll;
         LLGQ_arrays::atom_idx_y[atom] = 3 * atom_ll * n_coarse_ll + n_coarse_ll;
         LLGQ_arrays::atom_idx_z[atom] = 3 * atom_ll * n_coarse_ll + 2*n_coarse_ll;
      }
   }

   void calculate_random_fields(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse) {
      
      if (n_coarse < 0) {
         n_coarse = (n_fine > 0) ? ((n_fine - 1) / M + 1) : 0;
      }

      // Check if realizations is valid
      if (realizations <= 0 || n_fine <= 0 || M <= 0 || n_coarse <= 0) {
         std::cerr << "Error: Invalid parameters for quantum noise generation." << std::endl;
         std::cerr << "realizations=" << realizations << ", n_fine=" << n_fine 
                   << ", M=" << M << ", n_coarse=" << n_coarse << std::endl;
         return;
      }

      // Additional memory check with detailed output using 64-bit arithmetic
      const long long total_elements = static_cast<long long>(realizations) * static_cast<long long>(n_coarse);
      const long long memory_bytes = total_elements * sizeof(double);
      std::cout << "Memory allocation request: " << memory_bytes / (1024*1024) << " MB for " 
                << total_elements << " elements" << std::endl;
      
      // Check for potential overflow before allocation
      if (total_elements > static_cast<long long>(std::vector<double>().max_size())) {
         std::cerr << "Error: Requested vector size (" << total_elements 
                   << ") exceeds maximum allowed size (" << std::vector<double>().max_size() << ")" << std::endl;
         err::vexit();
      }

      // Check if FFTW is available
      #ifdef FFT
      if (n_fine <= 0) return; // Do not generate noise if there are no steps

      const double S0 = sim::internal::mp[0].S0.get();
      const double inv_sqrt_S0 = (S0 > 0) ? 1.0 / std::sqrt(S0) : 1.0;

      // Calculate coarse time step
      const double dt_coarse = dt_fine * M;

      double* __restrict in = (double*)fftw_malloc(sizeof(double) * n_coarse);
      fftw_complex* __restrict out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n_coarse/2 + 1));
      double* __restrict result = (double*)fftw_malloc(sizeof(double) * n_coarse);

      fftw_plan forward = fftw_plan_dft_r2c_1d(n_coarse, in, out, FFTW_MEASURE);
      fftw_plan backward = fftw_plan_dft_c2r_1d(n_coarse, out, result, FFTW_MEASURE);

      std::cout << "FFTW plans created for coarse grid." << std::endl;
      std::cout << "Total number fine time steps: " << n_fine << std::endl;
      std::cout << "Decimation factor M: " << M << std::endl;
      std::cout << "Number of coarse time steps: " << n_coarse << std::endl;
      std::cout << "FFT speedup factor: " << static_cast<double>(n_fine) / n_coarse << "x" << std::endl;
 
      const double norm_factor = 1.0 / n_coarse;

      // White noise statistics use fine time step to preserve variance
      static thread_local std::random_device rd;
      static thread_local std::mt19937 gen(rd());
      std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt_fine));

      try {
         // Use 64-bit arithmetic to avoid overflow
         const size_t safe_size = static_cast<size_t>(total_elements);
         LLGQ_arrays::coarse_noise_field.resize(safe_size);
         std::cout << "Successfully allocated noise field vector with size: " << safe_size 
                   << " (" << (safe_size*sizeof(double))/(1024*1024) << " MB)" << std::endl;
      } catch (const std::length_error& e) {
         std::cerr << "std::length_error during resize: " << e.what() << std::endl;
         std::cerr << "Requested size: " << total_elements << std::endl;
         std::cerr << "Max size: " << LLGQ_arrays::coarse_noise_field.max_size() << std::endl;
         err::vexit();
      } catch (const std::bad_alloc& e) {
         std::cerr << "std::bad_alloc during resize: " << e.what() << std::endl;
         err::vexit();
      }

      // Progress bar setup
      const int bar_width = 50;
      int last_printed_percent = -1;
      std::cout << "Generating PSD-based quantum noise fields on coarse grid..." << std::endl;

      // Precompute PSD on coarse grid frequency space
      std::vector<double> sqrt_PSD_coarse(n_coarse/2 + 1);
      double df_coarse = 1.0 / (n_coarse * dt_coarse);
      for (int i = 0; i <= n_coarse/2; ++i) {
         double omega = 2.0 * M_PI * i * df_coarse;
         sqrt_PSD_coarse[i] = std::sqrt(PSD(omega, T));
      }

      std::cout << "Starting noise generation for " << realizations << " realizations..." << std::endl;
      for (int r = 0; r < realizations; ++r) {
         // Generate white noise on the coarse grid
         for (int i = 0; i < n_coarse; ++i) {
            in[i] = dist(gen);
         }

         // Forward FFT
         fftw_execute(forward);

         // Apply PSD on coarse grid
         for (int i = 0; i <= n_coarse/2; i++) {
            const double magnitude = sqrt_PSD_coarse[i];
            out[i][0] *= magnitude;
            out[i][1] *= magnitude;
         }

         // Inverse FFT
         fftw_execute(backward);

         // Scale to preserve fine-time-step variance and store coarse noise
         const double scale = std::sqrt(dt_fine / dt_coarse);
         for (int j = 0; j < n_coarse; ++j) {
            // Use 64-bit arithmetic to prevent array index overflow
            const size_t index = static_cast<size_t>(j) + static_cast<size_t>(r) * static_cast<size_t>(n_coarse);
            LLGQ_arrays::coarse_noise_field[index] = result[j] * norm_factor * inv_sqrt_S0 * scale;
         }

         // Progress bar update - show progress every 5%
         int current_percent = static_cast<int>((r + 1) * 100.0 / realizations);
         if (current_percent >= last_printed_percent + 5 || r == realizations - 1) {
            float progress = static_cast<float>(r + 1) / realizations;
            int pos = static_cast<int>(bar_width * progress);

            std::cout << "\r[";
            for (int i = 0; i < bar_width; ++i) {
               if (i < pos) std::cout << "=";
               else if (i == pos) std::cout << ">";
               else std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << current_percent << "% (" << (r + 1) << "/" << realizations << ")";
            std::cout.flush();

            last_printed_percent = current_percent;
         }
      }

      // Cleanup FFTW resources after all realizations are complete
      fftw_destroy_plan(forward);
      fftw_destroy_plan(backward);
      fftw_free(in);
      fftw_free(out);
      fftw_free(result);

      std::cout << std::endl; // New line after progress bar

      #else
         std::cerr << "Error - quantum thermostat requires the FFTW library to function. Please recompile with the FFT library linked" << std::endl;
         err::vexit();
      #endif
   }

   double PSD(const double& omega, const double& T) {
      const double A = sim::internal::mp[0].A.get();
      const double Gamma = sim::internal::mp[0].Gamma.get();
      const double omega0 = sim::internal::mp[0].omega0.get();

      double x = (T > 1e-12) ? omega / (2 * T) : omega;  // Avoid division by zero
      double lorentzian_denom = (omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega;
      if (lorentzian_denom < 1e-12) lorentzian_denom = 1e-12; // Avoid division by zero
      double lorentzian = A * Gamma * omega / lorentzian_denom;
      double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero


      switch (sim::noise_type) {
         case 0: // Classical Noise
         return 2*T* A * Gamma / lorentzian_denom;
         case 1: // Quantum Noise
         if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
         else return coth * lorentzian;
         case 2: // Semiquantum Noise
         if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
         else return (coth-1) * lorentzian;

         default:
         std::cout << "Default Noise: Classical Noise" << std::endl;
         return 1.0 / x * lorentzian;
      }

   }

   void precompute_sqrt_PSD(int n, double dt, double T) {
      if (n <= 0) return;
      LLGQ_arrays::sqrt_PSD_buffer.resize(n/2 + 1);
      double df = 1.0 / (n * dt);
      for (int i = 0; i <= n/2; ++i) {
         double omega = 2.0 * M_PI * i * df;
         LLGQ_arrays::sqrt_PSD_buffer[i] = std::sqrt(PSD(omega, T));
      }
   }

   // Function to estimate cutoff frequency from noise_interpolated.cpp
   double estimate_cutoff_omega_cdf(double T, double target_frac) { 
       const double omega0 = sim::internal::mp[0].omega0.get();
       if (omega0 <= 0) return 1.0;
       double omega_max = 10.0 * omega0;
       int steps = 50000;
       std::vector<double> psd_vals(steps + 1);
       double domega = omega_max / steps;
       for (int i = 0; i <= steps; ++i) {
           double omega = i * domega;
           psd_vals[i] = PSD(omega, T);
       }
       double total_area = 0.0;
       for (int i = 0; i < steps; ++i) {
           total_area += 0.5 * (psd_vals[i] + psd_vals[i + 1]) * domega;
       }
       if (total_area <= 1e-12) return omega_max;
       double cum_area = 0.0;
       for (int i = 0; i <= steps; ++i) {
           if (i > 0) cum_area += 0.5 * (psd_vals[i - 1] + psd_vals[i]) * domega;
           if (cum_area >= target_frac * total_area) {
               return i * domega;
           }
       }
       return omega_max;
   }    

   

   double get_noise(const std::vector<double>& coarse_noise, double fine_step_idx, int M, size_t atom_idx) {
          double coarse_idx_float = fine_step_idx / M;
          size_t j = static_cast<size_t>(coarse_idx_float);
          double frac = coarse_idx_float - j;

          // Linear interpolation with safe index calculation
          const size_t index1 = j + atom_idx;
          const size_t index2 = j + atom_idx + 1;
          return coarse_noise[index1] * (1.0 - frac) + coarse_noise[index2] * frac;
    }

} // end of sim namespace
