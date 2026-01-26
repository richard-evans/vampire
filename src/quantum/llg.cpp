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
#include <vector>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "quantum.hpp"
#include "material.hpp"
#include "vio.hpp"

#include "internal.hpp"

#ifdef CUDA
#include "llg-HO-cuda.hpp"
#endif

namespace quantum{

   namespace internal{

      //---------------------------------------------------------------------------
      // Spin Dynamics Function
      //---------------------------------------------------------------------------
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

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // LLG Wrapper Function
   //---------------------------------------------------------------------------
   void llg_step(){

      using namespace internal;

      if(!initialised){
         // error
         zlog << zTs()  << "Programmer error: quantum noise is not initialised, exiting" << std::endl;
         std::cerr << "Programmer error: quantum noise is not initialised, exiting" << std::endl;
         err::vexit();
      }

      // Local variables
      const int num_atoms = atoms::num_atoms;
      std::vector<double> H(3);

      // Store initial state
      for (int atom = 0; atom < num_atoms; atom++) {
         y_in_storage[atom][0] = atoms::x_spin_array[atom];
         y_in_storage[atom][1] = atoms::y_spin_array[atom];
         y_in_storage[atom][2] = atoms::z_spin_array[atom];
         y_in_storage[atom][3] = x_v_array[atom];
         y_in_storage[atom][4] = y_v_array[atom];
         y_in_storage[atom][5] = z_v_array[atom];
         y_in_storage[atom][6] = x_w_array[atom];
         y_in_storage[atom][7] = y_w_array[atom];
         y_in_storage[atom][8] = z_w_array[atom];
      }


      // Calculate fields
      // Note: These functions are in sim namespace, but we need to access them.
      // They are likely not exported in a header we included?
      // sim::calculate_spin_fields and sim::calculate_external_fields are usually in sim.hpp or similar.
      // Let's assume they are available via sim.hpp.
      sim::calculate_spin_fields(0, num_atoms);
      sim::calculate_external_fields(0, num_atoms);

      const double dt = mp::dt;
      const double half_dt = 0.5 * dt;
      const double dt_over_6 = dt / 6.0;

      // K1 Step
      for (int atom = 0; atom < num_atoms; ++atom) {
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + quantum::get_field(atom, 0, 0.0);
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + quantum::get_field(atom, 1, 0.0);
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + quantum::get_field(atom, 2, 0.0);

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
      sim::calculate_spin_fields(0, num_atoms);

      // Call LLG method based on integration type
      #ifdef MPICF
         // MPI parallel version
         if(internal::llg_method == internal::llg_fft){
            internal::llg_FFT_mpi();
         }
         else if(internal::llg_method == internal::llg_ho){
            internal::llg_HO_mpi();
         }
      #else
         // Serial version
         #ifdef CUDA
            // CUDA accelerated version
            if(internal::llg_method == internal::llg_ho){
               internal::cuda_ho::llg_HO_step();
            }
            else if(internal::llg_method == internal::llg_fft){
               internal::llg_FFT();
            }
         #else
            // CPU serial version
            if(internal::llg_method == internal::llg_fft){
               internal::llg_FFT();
            }
            else if(internal::llg_method == internal::llg_ho){
               internal::llg_HO();
            }
         #endif
      #endif

      return;
   }

} // end of quantum namespace
