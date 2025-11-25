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
#include <cstring>

// Library for FFT
#ifdef FFT
#include <fftw3.h>
#endif

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"

#include "internal.hpp"

namespace quantum{
   namespace internal{

      //---------------------------------------------------------------------------
      // PSD function
      //---------------------------------------------------------------------------
      double PSD(const double& omega, const double& T) {
         // TODO: Support multiple materials? For now use material 0
         if(mp.empty()) return 0.0;

         const double A = mp[0].A.get();
         const double Gamma = mp[0].Gamma.get();
         const double omega0 = mp[0].omega0.get();

         double x = (T > 1e-12) ? omega / (2 * T) : omega;  // Avoid division by zero
         double lorentzian_denom = (omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega;
         if (lorentzian_denom < 1e-12) lorentzian_denom = 1e-12; // Avoid division by zero
         double lorentzian = A * Gamma * omega / lorentzian_denom;
         double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero

         switch (noise_type) {
            case 0: // Classical Noise
               return 2*T* A * Gamma / lorentzian_denom;
            case 1: // Quantum Noise
               if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
               else return coth * lorentzian;
            case 2: // Semiquantum Noise
               if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
               else return (coth-1) * lorentzian;
            default:
               return 1.0 / x * lorentzian;
         }
      }

      
      //---------------------------------------------------------------------------
      // Assign unique indices
      //---------------------------------------------------------------------------
      void assign_unique_indices(int n_coarse) {
            const int num_atoms = atoms::num_atoms;

            std::cout << "Assigning indices for " << num_atoms << " atoms with " << n_coarse << " coarse steps." << std::endl;

            atom_idx_x.resize(num_atoms);
            atom_idx_y.resize(num_atoms);
            atom_idx_z.resize(num_atoms);

            for (int atom = 0; atom < num_atoms; atom++) {
                // Use 64-bit arithmetic to prevent overflow
                const size_t atom_ll = static_cast<size_t>(atom);
                const size_t n_coarse_ll = static_cast<size_t>(n_coarse);
                atom_idx_x[atom] = 3 * atom_ll * n_coarse_ll;
                atom_idx_y[atom] = 3 * atom_ll * n_coarse_ll + n_coarse_ll;
                atom_idx_z[atom] = 3 * atom_ll * n_coarse_ll + 2*n_coarse_ll;
            }
      }

      //---------------------------------------------------------------------------
      // Windowed Random Field Calculation
      //---------------------------------------------------------------------------
      void calculate_noise(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse_total, std::vector<double>& noise_field) {
         #ifdef FFT
         if (window_size % 6 != 0) {
             std::cerr << "Error: Window size must be divisible by 6." << std::endl;
             return;
         }
         int segment_size = window_size / 6;

         if(mp.empty()) return;
         const double S0 = mp[0].S0.get();
         const double inv_sqrt_S0 = (S0 > 0) ? 1.0 / std::sqrt(S0) : 1.0;
         const double dt_coarse = dt_fine * M;
         const double norm_factor = 1.0 / window_size;
         const double scale = std::sqrt(dt_fine / dt_coarse);

         // FFTW Setup for the window size
         double* in = (double*)fftw_malloc(sizeof(double) * window_size);
         fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (window_size / 2 + 1));
         double* result = (double*)fftw_malloc(sizeof(double) * window_size);

         fftw_plan forward = fftw_plan_dft_r2c_1d(window_size, in, out, FFTW_MEASURE);
         fftw_plan backward = fftw_plan_dft_c2r_1d(window_size, out, result, FFTW_MEASURE);

         // Random number generator
         static thread_local std::random_device rd;
         static thread_local std::mt19937 gen(rd());
         std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt_fine));

         // Precompute PSD for the window
         std::vector<double> sqrt_PSD_window(window_size / 2 + 1);
         double df_window = 1.0 / (window_size * dt_coarse);
         for (int i = 0; i <= window_size / 2; ++i) {
             double omega = 2.0 * M_PI * i * df_window;
             sqrt_PSD_window[i] = std::sqrt(PSD(omega, T));
         }

         // Buffer to hold white noise
         std::vector<double> white_noise_buffer(window_size);

         // Resize output vector
         // Use 64-bit arithmetic for size calculation
         size_t total_size = static_cast<size_t>(n_coarse_total) * static_cast<size_t>(realizations);
         try {
            noise_field.resize(total_size);
         } catch (const std::bad_alloc& e) {
            std::cerr << "Error: Failed to allocate memory for noise field (" << total_size * sizeof(double) / (1024*1024) << " MB)" << std::endl;
            err::vexit();
         }

         std::cout << "Generating windowed noise..." << std::endl;

         for (int r = 0; r < realizations; ++r) {
             int generated_samples = 0;
             bool first_run = true;
             size_t realization_offset = static_cast<size_t>(r) * static_cast<size_t>(n_coarse_total);

             while (generated_samples < n_coarse_total) {

                 if (first_run) {
                     // Fill entire buffer with new random numbers
                     for (int i = 0; i < window_size; ++i) {
                         white_noise_buffer[i] = dist(gen);
                     }
                 } else {
                     // Shift buffer: Move last 2 segments (5 and 6) to front (1 and 2)
                     std::memmove(white_noise_buffer.data(),
                                  white_noise_buffer.data() + 4 * segment_size,
                                  2 * segment_size * sizeof(double));

                     // Fill the rest (segments 3, 4, 5, 6) with new random numbers
                     for (int i = 2 * segment_size; i < window_size; ++i) {
                         white_noise_buffer[i] = dist(gen);
                     }
                 }

                 // Copy to FFT input
                 for(int i=0; i<window_size; ++i) in[i] = white_noise_buffer[i];

                 // FFT
                 fftw_execute(forward);

                 // Apply PSD
                 for (int i = 0; i <= window_size / 2; i++) {
                     const double magnitude = sqrt_PSD_window[i];
                     out[i][0] *= magnitude;
                     out[i][1] *= magnitude;
                 }

                 // Inverse FFT
                 fftw_execute(backward);

                 // Process Output
                 if (first_run) {
                     // Take segments 1-5
                     int count = 5 * segment_size;
                     for (int j = 0; j < count; ++j) {
                         if (generated_samples < n_coarse_total) {
                              noise_field[realization_offset + generated_samples] = result[j] * norm_factor * inv_sqrt_S0 * scale;
                              generated_samples++;
                         }
                     }
                     first_run = false;
                 } else {
                     // Take segments 2-5
                     int start_idx = 1 * segment_size; // Start of segment 2
                     int end_idx = 5 * segment_size;   // End of segment 5 (exclusive)

                     for (int j = start_idx; j < end_idx; ++j) {
                         if (generated_samples < n_coarse_total) {
                             noise_field[realization_offset + generated_samples] = result[j] * norm_factor * inv_sqrt_S0 * scale;
                             generated_samples++;
                         }
                     }
                 }
             }
         }

         fftw_destroy_plan(forward);
         fftw_destroy_plan(backward);
         fftw_free(in);
         fftw_free(out);
         fftw_free(result);
         #else
         std::cerr << "Error: FFTW library required for quantum noise." << std::endl;
         err::vexit();
         #endif
      }

      //---------------------------------------------------------------------------
      // Get noise value with interpolation
      //---------------------------------------------------------------------------
      double get_noise(const std::vector<double>& coarse_noise, double fine_step_idx, int M, size_t atom_idx) {
          double coarse_idx_float = fine_step_idx / M;
          size_t j = static_cast<size_t>(coarse_idx_float);
          double frac = coarse_idx_float - j;

          // Linear interpolation with safe index calculation
          const size_t index1 = j + atom_idx;
          const size_t index2 = j + atom_idx + 1;

          // Boundary check (optional, but good for safety)
          if (index2 >= coarse_noise.size()) return coarse_noise[index1];

          return coarse_noise[index1] * (1.0 - frac) + coarse_noise[index2] * frac;
      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Get field wrapper
   //---------------------------------------------------------------------------
   double get_field(int atom, int component, double time_offset) {
       using namespace internal;
       size_t idx = 0;
       if(component == 0) idx = atom_idx_x[atom];
       else if(component == 1) idx = atom_idx_y[atom];
       else idx = atom_idx_z[atom];

       return get_noise(coarse_noise_field, noise_index + time_offset, M_decimation, idx);
   }

   //---------------------------------------------------------------------------
   // Increment time
   //---------------------------------------------------------------------------
   void increment_time() {
       internal::noise_index += 1.0;
   }

} // end of quantum namespace
