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
#include "vio.hpp"
#include "quantum.hpp"
#include "random.hpp"
#include "sim.hpp"

#include "internal.hpp"

namespace quantum{
   namespace internal{

      //---------------------------------------------------------------------------
      // PSD function
      //---------------------------------------------------------------------------
      double PSD(const double omega, const double T, const int material) {

         const double A      = material_A_array[material];
         const double Gamma  = material_gamma_array[material];
         const double omega0 = material_omega0_array[material];

         const double x = (T > 1e-12) ? omega / (2.0 * T) : omega;  // Avoid division by zero

         double lorentzian_denom = (omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega;
         if (lorentzian_denom < 1e-12) lorentzian_denom = 1e-12; // Avoid division by zero

         const double lorentzian = A * Gamma * omega / lorentzian_denom;
         const double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero (Joe Barker's idea)

         switch (internal::noise_type) {

            // Classical Noise
            case internal::classical:
               return 2.0 * T * A * Gamma / lorentzian_denom;
               break;

            // Quantum Noise
            case internal::quantum_zero:
               if(omega > 0.0) return coth * lorentzian;
               else return 2.0 * T * A * Gamma / (omega0 * omega0  * omega0 * omega0);
               break;

            // Quantum no-zero Noise
            case internal::quantum_no_zero:
               if(omega > 0.0) return (coth - 1.0) * lorentzian;
               else return 2.0 * T * A * Gamma / (omega0 * omega0  * omega0 * omega0);
               break;

            // Classical noise
            default:
               zlog << zTs()  << "Programmer error: quantum noise type is set to unknown value " << internal::noise_type << " which is not known" << std::endl;
               std::cerr << "Programmer error: quantum noise type is set to unknown value " << internal::noise_type << " which is not known" << std::endl;
               err::vexit();

         }
      }
      
   }


      //---------------------------------------------------------------------------
      // Assign unique indices to map noise for each atom and spatial dimension
      //---------------------------------------------------------------------------
      void assign_unique_indices(int n_coarse, int num_atoms_local) {
            std::cout << "Assigning indices for " << num_atoms_local << " local atoms with " << n_coarse << " coarse steps." << std::endl;

            atom_idx_x.resize(num_atoms_local);
            atom_idx_y.resize(num_atoms_local);
            atom_idx_z.resize(num_atoms_local);

            for (int atom = 0; atom < num_atoms_local; atom++) {
                // Use 64-bit arithmetic to prevent overflow
                const size_t atom_ll = static_cast<size_t>(atom);
                const size_t n_coarse_ll = static_cast<size_t>(n_coarse);
                atom_idx_x[atom] = 3 * atom_ll * n_coarse_ll;
                atom_idx_y[atom] = 3 * atom_ll * n_coarse_ll + n_coarse_ll;
                atom_idx_z[atom] = 3 * atom_ll * n_coarse_ll + 2*n_coarse_ll;
            }
      }

      //---------------------------------------------------------------------------
      // Non-windowed Random Field Calculation (Direct FFT on coarse grid)
      //---------------------------------------------------------------------------
      void calculate_noise(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse_total, std::vector<double>& noise_field) {
         #ifdef FFT
         
         if (n_fine <= 0) return; // Do not generate noise if there are no steps
         
         const double dt_coarse = dt_fine * M;

         // Allocate FFTW arrays for the full coarse grid
         double* __restrict in = (double*)fftw_malloc(sizeof(double) * n_coarse_total);
         fftw_complex* __restrict out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n_coarse_total/2 + 1));
         double* __restrict result = (double*)fftw_malloc(sizeof(double) * n_coarse_total);

         // Create FFTW plans
         fftw_plan forward = fftw_plan_dft_r2c_1d(n_coarse_total, in, out, FFTW_MEASURE);
         fftw_plan backward = fftw_plan_dft_c2r_1d(n_coarse_total, out, result, FFTW_MEASURE);

         std::cout << "FFTW plans created for coarse grid." << std::endl;
         std::cout << "Total number fine time steps: " << n_fine << std::endl;
         std::cout << "Decimation factor M: " << M << std::endl;
         std::cout << "Number of coarse time steps: " << n_coarse_total << std::endl;
         std::cout << "FFT speedup factor: " << static_cast<double>(n_fine) / n_coarse_total << "x" << std::endl;

         const double norm_factor = 1.0 / n_coarse_total;

         // White noise statistics use fine time step to preserve variance
         // Each MPI process generates independent noise for its local atoms
         static thread_local std::random_device rd;
         static thread_local std::mt19937 gen(rd());
         std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt_fine));

         // Resize output vector with 64-bit arithmetic to avoid overflow
         try {
            const size_t total_elements = static_cast<size_t>(realizations) * static_cast<size_t>(n_coarse_total);
            noise_field.resize(total_elements);
            std::cout << "Successfully allocated noise field vector with size: " << total_elements 
                     << " (" << (total_elements*sizeof(double))/(1024*1024) << " MB)" << std::endl;
         } catch (const std::length_error& e) {
            std::cerr << "std::length_error during resize: " << e.what() << std::endl;
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
         std::vector<double> sqrt_PSD_coarse(n_coarse_total/2 + 1);
         double df_coarse = 1.0 / (n_coarse_total * dt_coarse);
         for (int i = 0; i <= n_coarse_total/2; ++i) {
            double omega = 2.0 * M_PI * i * df_coarse;
            sqrt_PSD_coarse[i] = std::sqrt(PSD(omega, T, 0)); // Use material 0 for now
         }

         // Scale factor to preserve fine-time-step variance
         const double scale = std::sqrt(dt_fine / dt_coarse);

         std::cout << "Starting noise generation for " << realizations << " realizations..." << std::endl;
         for (int r = 0; r < realizations; ++r) {
            // Generate white noise on the coarse grid
            for (int i = 0; i < n_coarse_total; ++i) {
               in[i] = dist(gen);
            }

            // Forward FFT
            fftw_execute(forward);

            // Apply PSD on coarse grid
            for (int i = 0; i <= n_coarse_total/2; i++) {
               const double magnitude = sqrt_PSD_coarse[i];
               out[i][0] *= magnitude;
               out[i][1] *= magnitude;
            }

            // Inverse FFT
            fftw_execute(backward);

            // Store coarse noise with proper scaling
            for (int j = 0; j < n_coarse_total; ++j) {
               // Use 64-bit arithmetic to prevent array index overflow
               const size_t index = static_cast<size_t>(j) + static_cast<size_t>(r) * static_cast<size_t>(n_coarse_total);
               noise_field[index] = result[j] * norm_factor * scale;
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

   namespace internal {

      //---------------------------------------------------------------------------
      // Initialize noise structures (stub - to be implemented if needed)
      //---------------------------------------------------------------------------
      void init_noise_structures(int n_coarse, int realizations, double dt_fine, int M, double T) {
         // This function should initialize FFT noise structures
         // For now, just a stub to allow compilation
         return;
      }

   } // end of internal namespace

} // end of quantum namespace
