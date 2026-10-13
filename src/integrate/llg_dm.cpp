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
//   ported from src/jams/solvers/cuda_llg_dm.cu
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
//   Depondt-Mertens LLG integrator (llg-dm).
//
//   Deterministic predictor-corrector step:
//      k1 = dt omega(S_n);  S* = R(k1) S_n
//      k  = (k1 + dt omega(S*)) / 2;  S_n+1 = R(k) S_n
//
//   Thermal noise is applied by Strang splitting: a dt/2 stochastic rotation
//   before and after the deterministic step, each with a fresh random field.
//
//   Reference:
//      Ph. Depondt and F. G. Mertens, J. Phys.: Condens. Matter 21, 336005 (2009)
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

namespace integrate{
namespace internal{

//------------------------------------------------------------------------------
// Predictor stage: store S_n, k1 = dt omega(S_n), S* = R(k1) S_n
//------------------------------------------------------------------------------
void dm_predictor_stage(const int start_index, const int end_index){

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

      llg_rotate_initial_spin(atom, k1);

   }

   return;

}

//------------------------------------------------------------------------------
// Corrector stage: k = (k1 + dt omega(S*)) / 2, S_n+1 = R(k) S_n
//------------------------------------------------------------------------------
void dm_corrector_stage(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){

      double k2[3];
      llg_rotation_vector(atom, k2);

      const double k[3] = { 0.5 * (x_k1_array[atom] + k2[0]),
                            0.5 * (y_k1_array[atom] + k2[1]),
                            0.5 * (z_k1_array[atom] + k2[2]) };

      llg_rotate_initial_spin(atom, k);

   }

   return;

}

//------------------------------------------------------------------------------
// Function to integrate the system by one time step using Depondt-Mertens
//------------------------------------------------------------------------------
void llg_dm_step(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "integrate::internal::llg_dm_step has been called" << std::endl;}

   // Check for initialisation of integration arrays
   if(LLG_rotation_arrays::LLG_rotation_set == false) LLG_rotation_init();

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   // deterministic step (external fields calculated once per step)
   llg_rotation_stage(dm_predictor_stage, true);
   llg_rotation_stage(dm_corrector_stage, false);

   // noise half step
   llg_rotation_noise_half_step(0, atoms::num_atoms);

   return;

}

} // end of internal namespace
} // end of integrate namespace
