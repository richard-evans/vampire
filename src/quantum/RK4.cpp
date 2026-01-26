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
#include "errors.hpp"
#include "vio.hpp"

#include "internal.hpp"

namespace quantum{
   namespace internal{

      // Equation of motion for spin and Kernel auxiliary variables
      void spinDynamics(const double* y, const double* H, double* dydt, const int material) {

            // Get lorentzian parameters
            const double A      = material_A_array[material];
            const double Gamma  = material_gamma_array[material];
            const double omega0 = material_omega0_array[material];

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

      void LL_HO_method(const double* y, const double* H, double* dydt, const int material, const double dt) {
          // Get harmonic oscillator parameters
          const double A      = material_A_array[material];
          const double Gamma  = material_gamma_array[material];
          const double omega0 = material_omega0_array[material];
          const double S0     = material_S0_array[material];
          const double T      = sim::temperature * 1.34; // Convert from K to energy units

          // Calculate noise prefactor once
          const double noise_prefactor = std::sqrt(4.0 * Gamma * A * T / dt);

          // Generate thermal noise
          const double noise_x = noise_prefactor * mtrandom::gaussian();
          const double noise_y = noise_prefactor * mtrandom::gaussian();
          const double noise_z = noise_prefactor * mtrandom::gaussian();

          // y[0-2]: spin components (S)
          // y[3-5]: position components (q)
          // y[6-8]: momentum components (w/p)

          // dS/dt = S × (H + q)
          dydt[0] =  (y[1]*(H[2]+y[5]) - y[2]*(H[1]+y[4]));
          dydt[1] =  (y[2]*(H[0]+y[3]) - y[0]*(H[2]+y[5]));
          dydt[2] =  (y[0]*(H[1]+y[4]) - y[1]*(H[0]+y[3]));

          // dq/dt = p 
          dydt[3] = y[6];
          dydt[4] = y[7];
          dydt[5] = y[8];

          // dp/dt = -ω₀²q - Γw + AS + noise
          dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0] + noise_x;
          dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1] + noise_y;
          dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2] + noise_z;
      }

      // FFT version: same 9-component equations as HO, but H already includes thermal noise
      // Thermal noise is added to H externally before calling this function
      void LL_FFT_method(const double* y, const double* H, double* dydt, const int material) {
          // Get harmonic oscillator parameters
          const double A      = material_A_array[material];
          const double Gamma  = material_gamma_array[material];
          const double omega0 = material_omega0_array[material];

          // y[0-2]: spin components (S)
          // y[3-5]: position components (q - or V in Paper)
          // y[6-8]: momentum components (p - or W in Paper)

          // dS/dt = S × (H + q)
          dydt[0] = y[1]*(H[2] + y[5]) - y[2]*(H[1] + y[4]);
          dydt[1] = y[2]*(H[0] + y[3]) - y[0]*(H[2] + y[5]);
          dydt[2] = y[0]*(H[1] + y[4]) - y[1]*(H[0] + y[3]);

          // dq/dt = p
          dydt[3] = y[6];
          dydt[4] = y[7];
          dydt[5] = y[8];

          // dp/dt = -ω₀²q - Γp + AS (NO noise term - it's in H)
          dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0];
          dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1];
          dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2];
      }

   } // end of internal namespace
} // end of quantum namespace