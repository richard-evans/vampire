//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
//   This file is a derivative work of JAMS (https://github.com/stonerlab/jams),
//   ported from src/jams/solvers/cpu_llg_rkmk.cc, cuda_llg_rkmk2.cu and cuda_llg_rkmk4.cu
//   (version3 branch, commit 234ba11).
//   Original implementation by Dr Joseph Barker, University of Leeds.
//
//   The original JAMS code is distributed under the MIT License, reproduced
//   below in full as required:
//
//   MIT License
//
//   Copyright (c) 2009-2013 Joseph Barker (University of York, UK),
//                 2014-2018 Joseph Barker (Tohoku University, Japan),
//                 2018-2024 Joseph Barker and collaborators (University of Leeds, UK)
//
//   Permission is hereby granted, free of charge, to any person obtaining a copy
//   of this software and associated documentation files (the "Software"), to deal
//   in the Software without restriction, including without limitation the rights
//   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//   copies of the Software, and to permit persons to whom the Software is
//   furnished to do so, subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in all
//   copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.
//
//------------------------------------------------------------------------------
//
//   Runge-Kutta-Munthe-Kaas LLG integrators (llg-rkmk2 and llg-rkmk4).
//
//   Deterministic step with rotation vectors k_i = dt omega_i:
//
//   RKMK2 (explicit midpoint tableau)
//      k1 = dt omega(S0);  S1/2 = R(k1/2) S0
//      k  = dexp_inv(k1/2, dt omega(S1/2));  S1 = R(k) S0
//
//   RKMK4 (classical RK4 tableau)
//      k1 = dt omega(S0);                            S = R(k1/2) S0
//      k2 = dexp_inv(k1/2, dt omega(S));             S = R(k2/2) S0
//      k3 = dexp_inv(k2/2, dt omega(S));             S = R(k3) S0
//      k4 = dexp_inv(k3,   dt omega(S));
//      S1 = R((k1 + 2 k2 + 2 k3 + k4)/6) S0
//
//   Thermal noise is applied by Strang splitting: a dt/2 stochastic rotation
//   before and after the deterministic step, each with a fresh random field.
//
//   References:
//      H. Munthe-Kaas, BIT Numerical Mathematics 38, 92 (1998)
//      H. Munthe-Kaas, Applied Numerical Mathematics 29, 115 (1999)
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"

// integrate module headers
#include "integrate.hpp"
#include "internal.hpp"
#include "llg_rotation.hpp"
#include "rotations.hpp"

namespace integrate{
namespace internal{

//------------------------------------------------------------------------------
// First stage (RKMK2 and RKMK4): store S0, k1 = dt omega(S0), S = R(k1/2) S0
//------------------------------------------------------------------------------
void rkmk_stage_1(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      x_initial_spin_array[atom] = atoms::x_spin_array[atom];
      y_initial_spin_array[atom] = atoms::y_spin_array[atom];
      z_initial_spin_array[atom] = atoms::z_spin_array[atom];

      double k1[3];
      llg_rotation_vector(atom, k1);

      x_k1_array[atom] = k1[0];
      y_k1_array[atom] = k1[1];
      z_k1_array[atom] = k1[2];

      const double phi[3] = { 0.5 * k1[0], 0.5 * k1[1], 0.5 * k1[2] };
      llg_rotate_initial_spin(atom, phi);

   }

   return;

}

//------------------------------------------------------------------------------
// Final RKMK2 stage: k = dexp_inv(k1/2, dt omega(S1/2)), S1 = R(k) S0
//------------------------------------------------------------------------------
void rkmk2_stage_2(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      double v2[3];
      llg_rotation_vector(atom, v2);

      const double phi[3] = { 0.5 * x_k1_array[atom], 0.5 * y_k1_array[atom], 0.5 * z_k1_array[atom] };

      double k[3];
      rotations::dexp_inv_so3(phi, v2, k);

      llg_rotate_initial_spin(atom, k);

   }

   return;

}

//------------------------------------------------------------------------------
// Second RKMK4 stage: k2 = dexp_inv(k1/2, dt omega(S)), S = R(k2/2) S0
//------------------------------------------------------------------------------
void rkmk4_stage_2(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      double v2[3];
      llg_rotation_vector(atom, v2);

      const double phi[3] = { 0.5 * x_k1_array[atom], 0.5 * y_k1_array[atom], 0.5 * z_k1_array[atom] };

      double k2[3];
      rotations::dexp_inv_so3(phi, v2, k2);

      x_k2_array[atom] = k2[0];
      y_k2_array[atom] = k2[1];
      z_k2_array[atom] = k2[2];

      const double phi2[3] = { 0.5 * k2[0], 0.5 * k2[1], 0.5 * k2[2] };
      llg_rotate_initial_spin(atom, phi2);

   }

   return;

}

//------------------------------------------------------------------------------
// Third RKMK4 stage: k3 = dexp_inv(k2/2, dt omega(S)), S = R(k3) S0
//------------------------------------------------------------------------------
void rkmk4_stage_3(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      double v3[3];
      llg_rotation_vector(atom, v3);

      const double phi[3] = { 0.5 * x_k2_array[atom], 0.5 * y_k2_array[atom], 0.5 * z_k2_array[atom] };

      double k3[3];
      rotations::dexp_inv_so3(phi, v3, k3);

      x_k3_array[atom] = k3[0];
      y_k3_array[atom] = k3[1];
      z_k3_array[atom] = k3[2];

      llg_rotate_initial_spin(atom, k3);

   }

   return;

}

//------------------------------------------------------------------------------
// Final RKMK4 stage: k4 = dexp_inv(k3, dt omega(S)),
// S1 = R((k1 + 2 k2 + 2 k3 + k4)/6) S0
//------------------------------------------------------------------------------
void rkmk4_stage_4(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      double v4[3];
      llg_rotation_vector(atom, v4);

      const double phi[3] = { x_k3_array[atom], y_k3_array[atom], z_k3_array[atom] };

      double k4[3];
      rotations::dexp_inv_so3(phi, v4, k4);

      const double k[3] = { (x_k1_array[atom] + 2.0 * x_k2_array[atom] + 2.0 * x_k3_array[atom] + k4[0]) / 6.0,
                            (y_k1_array[atom] + 2.0 * y_k2_array[atom] + 2.0 * y_k3_array[atom] + k4[1]) / 6.0,
                            (z_k1_array[atom] + 2.0 * z_k2_array[atom] + 2.0 * z_k3_array[atom] + k4[2]) / 6.0 };

      llg_rotate_initial_spin(atom, k);

   }

   return;

}

//------------------------------------------------------------------------------
// Function to integrate the system by one time step using RKMK2
//------------------------------------------------------------------------------
void llg_rkmk2_step(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "integrate::internal::llg_rkmk2_step has been called" << std::endl;}

   // Check for initialisation of integration arrays
   if(LLG_rotation_arrays::LLG_rotation_set == false) LLG_rotation_init();

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   // deterministic step (external fields calculated once per step)
   llg_rotation_stage(rkmk_stage_1, true);
   llg_rotation_stage(rkmk2_stage_2, false);

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   return;

}

//------------------------------------------------------------------------------
// Function to integrate the system by one time step using RKMK4
//------------------------------------------------------------------------------
void llg_rkmk4_step(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "integrate::internal::llg_rkmk4_step has been called" << std::endl;}

   // Check for initialisation of integration arrays
   if(LLG_rotation_arrays::LLG_rotation_set == false) LLG_rotation_init();

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   // deterministic step (external fields calculated once per step)
   llg_rotation_stage(rkmk_stage_1, true);
   llg_rotation_stage(rkmk4_stage_2, false);
   llg_rotation_stage(rkmk4_stage_3, false);
   llg_rotation_stage(rkmk4_stage_4, false);

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   return;

}

} // end of internal namespace
} // end of integrate namespace
