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
//   ported from src/jams/solvers/cpu_llg_rkmk.cc
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
//   Shared storage, Strang-split noise substep and field evaluation for the
//   rotation based LLG integrators (llg-rkmk2, llg-rkmk4 and llg-dm).
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "gpu.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// integrate module headers
#include "integrate.hpp"
#include "internal.hpp"
#include "rotations.hpp"

namespace LLG_rotation_arrays{

   // Local arrays for rotation based LLG integration
   std::vector<double> x_initial_spin_array;
   std::vector<double> y_initial_spin_array;
   std::vector<double> z_initial_spin_array;

   std::vector<double> x_k1_array;
   std::vector<double> y_k1_array;
   std::vector<double> z_k1_array;

   std::vector<double> x_k2_array;
   std::vector<double> y_k2_array;
   std::vector<double> z_k2_array;

   std::vector<double> x_k3_array;
   std::vector<double> y_k3_array;
   std::vector<double> z_k3_array;

   std::vector<double> x_spin_storage_array;
   std::vector<double> y_spin_storage_array;
   std::vector<double> z_spin_storage_array;

   bool LLG_rotation_set = false; // Flag to define state of arrays (initialised/uninitialised)

}

namespace integrate{
namespace internal{

//------------------------------------------------------------------------------
// Function to check that the selected rotation based integrator is supported
// for this simulation and to log the method references
//------------------------------------------------------------------------------
void initialize_rotation_integrators(){

   const bool rkmk = integrate::integrator == integrator_t::llg_rkmk2 || integrate::integrator == integrator_t::llg_rkmk4;
   const bool dm   = integrate::integrator == integrator_t::llg_dm;

   if(!rkmk && !dm) return;

   // check for unsupported combinations
   std::string unsupported = "";
   if(gpu::acceleration)                                       unsupported = "GPU acceleration";
   if(micromagnetic::discretisation_type > 0)                  unsupported = "micromagnetic multiscale simulations";
   if(program::program == program::program_t::hamr_simulation) unsupported = "the HAMR program";

   if(unsupported != ""){
      std::ostringstream message;
      message << "Error - the selected rotation based LLG integrator (llg-rkmk2, llg-rkmk4 or llg-dm) is not yet supported with " << unsupported << ". Exiting.";
      err::fatal(message.str());
   }

   zlog << zTs() << "Rotation based LLG integrators ported from JAMS (https://github.com/stonerlab/jams), original implementation by Dr Joseph Barker, University of Leeds" << std::endl;
   if(rkmk){
      zlog << zTs() << "Runge-Kutta-Munthe-Kaas integrator: H. Munthe-Kaas, BIT 38, 92 (1998) and Appl. Numer. Math. 29, 115 (1999)" << std::endl;
   }
   if(dm){
      zlog << zTs() << "Depondt-Mertens integrator: Ph. Depondt and F. G. Mertens, J. Phys.: Condens. Matter 21, 336005 (2009)" << std::endl;
   }

   return;

}

//------------------------------------------------------------------------------
// Function to allocate arrays for rotation based integrators
//------------------------------------------------------------------------------
int LLG_rotation_init(){

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout << "integrate::internal::LLG_rotation_init has been called" << std::endl;}

   using namespace LLG_rotation_arrays;

   x_initial_spin_array.resize(atoms::num_atoms,0.0);
   y_initial_spin_array.resize(atoms::num_atoms,0.0);
   z_initial_spin_array.resize(atoms::num_atoms,0.0);

   x_k1_array.resize(atoms::num_atoms,0.0);
   y_k1_array.resize(atoms::num_atoms,0.0);
   z_k1_array.resize(atoms::num_atoms,0.0);

   x_k2_array.resize(atoms::num_atoms,0.0);
   y_k2_array.resize(atoms::num_atoms,0.0);
   z_k2_array.resize(atoms::num_atoms,0.0);

   x_k3_array.resize(atoms::num_atoms,0.0);
   y_k3_array.resize(atoms::num_atoms,0.0);
   z_k3_array.resize(atoms::num_atoms,0.0);

   x_spin_storage_array.resize(atoms::num_atoms,0.0);
   y_spin_storage_array.resize(atoms::num_atoms,0.0);
   z_spin_storage_array.resize(atoms::num_atoms,0.0);

   LLG_rotation_set = true;

   return EXIT_SUCCESS;

}

//------------------------------------------------------------------------------
// Function to apply a half time step stochastic rotation to spins in the
// range start_index - end_index (Strang splitting of the noise)
//
// The thermal field sigma (mp::material[].H_th_sigma) is set for the full
// time step dt; the noise substep has length dt/2 and so needs a field
// larger by sqrt(2) (sigma ~ 1/sqrt(dt)).
//------------------------------------------------------------------------------
void llg_rotation_noise_half_step(const int start_index, const int end_index){

   // generate new thermal fields
   sim::calculate_thermal_fields_only(start_index, end_index);

   const double noise_scale = std::sqrt(2.0);

   for(int atom = start_index; atom < end_index; atom++){

      const int imaterial = atoms::type_array[atom];
      const double gyro  = -mp::material[imaterial].one_oneplusalpha_sq;
      const double alpha = mp::material[imaterial].alpha;

      const double s[3] = { atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom] };
      const double h_noise[3] = { noise_scale * atoms::x_thermal_field_array[atom],
                                  noise_scale * atoms::y_thermal_field_array[atom],
                                  noise_scale * atoms::z_thermal_field_array[atom] };

      double s_out[3];
      rotations::noise_step_rodrigues(s, h_noise, gyro, alpha, mp::half_dt, s_out);

      atoms::x_spin_array[atom] = s_out[0];
      atoms::y_spin_array[atom] = s_out[1];
      atoms::z_spin_array[atom] = s_out[2];

   }

   return;

}

//------------------------------------------------------------------------------
// Function to copy new spins from storage to spin array
//------------------------------------------------------------------------------
void llg_rotation_copy_storage_to_spins(const int start_index, const int end_index){

   using namespace LLG_rotation_arrays;

   for(int atom = start_index; atom < end_index; atom++){
      atoms::x_spin_array[atom] = x_spin_storage_array[atom];
      atoms::y_spin_array[atom] = y_spin_storage_array[atom];
      atoms::z_spin_array[atom] = z_spin_storage_array[atom];
   }

   return;

}

//------------------------------------------------------------------------------
// Function to evaluate fields and apply one integrator stage to all atoms
//
// The stage function reads the current spins and deterministic fields and
// writes new spins to the storage arrays. External (non-thermal) fields are
// optionally recalculated, typically once at the start of each step.
//------------------------------------------------------------------------------
void llg_rotation_stage(void (*stage)(const int, const int), const bool calculate_external_fields){

   const int num_atoms = atoms::num_atoms;

   sim::calculate_spin_fields(0, num_atoms);
   if(calculate_external_fields) sim::calculate_deterministic_external_fields(0, num_atoms);

   stage(0, num_atoms);

   llg_rotation_copy_storage_to_spins(0, num_atoms);

   return;

}

#ifdef MPICF
//------------------------------------------------------------------------------
// Parallel version of llg_rotation_stage. Core atoms are updated while the
// halo swap is in progress, then boundary atoms once it has completed. New
// spins are only copied back once all fields have been evaluated.
//------------------------------------------------------------------------------
void llg_rotation_stage_mpi(void (*stage)(const int, const int), const bool calculate_external_fields){

   // Local variables for core / boundary integration
   const int pre_comm_si = 0;
   const int pre_comm_ei = vmpi::num_core_atoms;
   const int post_comm_si = vmpi::num_core_atoms;
   const int post_comm_ei = vmpi::num_core_atoms+vmpi::num_bdry_atoms;

   // Initiate halo swap
   vmpi::mpi_init_halo_swap();

   // Calculate fields and stage (core)
   sim::calculate_spin_fields(pre_comm_si, pre_comm_ei);
   if(calculate_external_fields) sim::calculate_deterministic_external_fields(pre_comm_si, pre_comm_ei);
   stage(pre_comm_si, pre_comm_ei);

   // Complete halo swap
   vmpi::mpi_complete_halo_swap();

   // Calculate fields and stage (boundary)
   sim::calculate_spin_fields(post_comm_si, post_comm_ei);
   if(calculate_external_fields) sim::calculate_deterministic_external_fields(post_comm_si, post_comm_ei);
   stage(post_comm_si, post_comm_ei);

   // Copy new spins to spin array (all)
   llg_rotation_copy_storage_to_spins(pre_comm_si, post_comm_ei);

   return;

}
#endif

} // end of internal namespace
} // end of integrate namespace
