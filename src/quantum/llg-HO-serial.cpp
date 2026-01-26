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
#include "random.hpp"

#include "internal.hpp"

namespace quantum{

   namespace internal{


   

   //---------------------------------------------------------------------------
   // LLG HO Step Function
   //---------------------------------------------------------------------------
   void llg_HO(){

      using namespace quantum;
      using namespace internal;

      // Local variables
      const int num_atoms = atoms::num_atoms;
      std::vector<double> H(3);
      /// @brief Position vector for the quantum harmonic oscillator
      /// @details Three-dimensional vector storing the position coordinates (q_x, q_y, q_z)
      ///          of the quantum harmonic oscillator in the LLG integration scheme
      // Local q and p arrays are already available as:
      // q_x_array, q_y_array, q_z_array (position)
      // p_x_array, p_y_array, p_z_array (momentum)
      /// @details Three-dimensional vector storing the momentum coordinates (p_x, p_y, p_z)
      ///          of the quantum harmonic oscillator in the LLG integration scheme
      std::vector<double> p_vec(3);

      // Local time step variables
      const double dt = mp::dt;
      const double half_dt = 0.5 * dt;
      const double dt_over_6 = dt / 6.0;

      // Store initial state (S, q, p)
      for (int atom = 0; atom < num_atoms; atom++) {
         // Spin components (0-2)
         y_in_storage[atom][0] = atoms::x_spin_array[atom];
         y_in_storage[atom][1] = atoms::y_spin_array[atom];
         y_in_storage[atom][2] = atoms::z_spin_array[atom];
         // Position components (3-5)
         y_in_storage[atom][3] = q_x_array[atom];
         y_in_storage[atom][4] = q_y_array[atom];
         y_in_storage[atom][5] = q_z_array[atom];
         // Momentum components (6-8)
         y_in_storage[atom][6] = p_x_array[atom];
         y_in_storage[atom][7] = p_y_array[atom];
         y_in_storage[atom][8] = p_z_array[atom];
      }

      // Calculate fields
      sim::calculate_spin_fields(0, num_atoms);
      sim::calculate_external_fields(0, num_atoms);

      // Now update both spins and lattice dynamics using RK4
      
      // K1 Step
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];

         // Calculate K1 for full state (S, q, p)
         LL_HO_method(y_in_storage[atom].data(), H.data(), k1_storage[atom].data(), imaterial, dt);

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

      // K2 Step
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];

         // Calculate K2
         LL_HO_method(y_pred_storage[atom].data(), H.data(), k2_storage[atom].data(), imaterial, dt);

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
      sim::calculate_spin_fields(0, num_atoms);

      // K3 Step
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];

         // Calculate K3
         LL_HO_method(y_pred_storage[atom].data(), H.data(), k3_storage[atom].data(), imaterial, dt);

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
      sim::calculate_spin_fields(0, num_atoms);

      //K4 step
      for (int atom = 0; atom < num_atoms; ++atom) {
         const int imaterial = atoms::type_array[atom];
         H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom];
         H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom];
         H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom];

         // Calculate K4
         LL_HO_method(y_pred_storage[atom].data(), H.data(), k4_storage[atom].data(), imaterial, dt);

         // Final update for all 9 components (S, q, p)
         for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + dt_over_6 * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
         }

         // Normalize spin components only (0-2)
         double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
         const double inv_mag = 1.0 / S_magnitude;
         y_pred_storage[atom][0] *= inv_mag;
         y_pred_storage[atom][1] *= inv_mag;
         y_pred_storage[atom][2] *= inv_mag;

         // Update spin arrays
         atoms::x_spin_array[atom] = y_pred_storage[atom][0];
         atoms::y_spin_array[atom] = y_pred_storage[atom][1];
         atoms::z_spin_array[atom] = y_pred_storage[atom][2];
         
         // Update lattice position arrays (q)
         q_x_array[atom] = y_pred_storage[atom][3];
         q_y_array[atom] = y_pred_storage[atom][4];
         q_z_array[atom] = y_pred_storage[atom][5];
         
         // Update lattice momentum arrays (p)
         p_x_array[atom] = y_pred_storage[atom][6];
         p_y_array[atom] = y_pred_storage[atom][7];
         p_z_array[atom] = y_pred_storage[atom][8];
      }

      return;
   }

} // end of internal namespace

} // end of quantum namespace
