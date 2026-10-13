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
//   Parallel (MPI) Depondt-Mertens LLG integrator (llg-dm).
//   Stages are shared with the serial version; see llg_rotation.cpp for the
//   core / boundary field evaluation around the halo swap.
//
//------------------------------------------------------------------------------
//
#ifdef MPICF

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "vmpi.hpp"

// integrate module headers
#include "integrate.hpp"
#include "internal.hpp"

namespace integrate{
namespace internal{

//------------------------------------------------------------------------------
// Function to integrate the system by one time step using Depondt-Mertens
// (parallel)
//------------------------------------------------------------------------------
void llg_dm_mpi_step(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "integrate::internal::llg_dm_mpi_step has been called" << std::endl;}

   // Check for initialisation of integration arrays
   if(LLG_rotation_arrays::LLG_rotation_set == false) LLG_rotation_init();

   // number of local atoms (core and boundary)
   const int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

   // noise half step
   llg_rotation_noise_half_step(0, num_local_atoms);

   // deterministic step (external fields calculated once per step)
   llg_rotation_stage_mpi(dm_predictor_stage, true);
   llg_rotation_stage_mpi(dm_corrector_stage, false);

   // noise half step
   llg_rotation_noise_half_step(0, num_local_atoms);

   // Wait for other processors
   vmpi::barrier();

   return;

}

} // end of internal namespace
} // end of integrate namespace
#endif
