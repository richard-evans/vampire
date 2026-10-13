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

#ifndef INTEGRATE_LLG_ROTATION_H_
#define INTEGRATE_LLG_ROTATION_H_

//------------------------------------------------------------------------------
// Inline per-atom helpers shared by the rotation based LLG integrators
// (llg-rkmk2, llg-rkmk4 and llg-dm), kept in a header so they are inlined
// into the stage loops
//------------------------------------------------------------------------------

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"

// integrate module headers
#include "internal.hpp"
#include "rotations.hpp"

namespace integrate{
namespace internal{

//------------------------------------------------------------------------------
// Rotation vector dt omega(S, H) for an atom using current spin and
// deterministic fields (spin + external, excluding thermal)
//------------------------------------------------------------------------------
inline void llg_rotation_vector(const int atom, double k[3]){

   const int imaterial = atoms::type_array[atom];
   const double gyro  = -mp::material[imaterial].one_oneplusalpha_sq;
   const double alpha = mp::material[imaterial].alpha;

   const double S[3] = { atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom] };
   const double H[3] = { atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                         atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                         atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] };

   double omega[3];
   rotations::omega_llg(S, H, gyro, alpha, omega);

   k[0] = mp::dt * omega[0];
   k[1] = mp::dt * omega[1];
   k[2] = mp::dt * omega[2];

}

//------------------------------------------------------------------------------
// Rotate initial spin S0 of atom by phi and write to storage array
//------------------------------------------------------------------------------
inline void llg_rotate_initial_spin(const int atom, const double phi[3]){

   using namespace LLG_rotation_arrays;

   const double S0[3] = { x_initial_spin_array[atom], y_initial_spin_array[atom], z_initial_spin_array[atom] };

   double S_new[3];
   rotations::rodrigues_rotate(phi, S0, S_new);

   x_spin_storage_array[atom] = S_new[0];
   y_spin_storage_array[atom] = S_new[1];
   z_spin_storage_array[atom] = S_new[2];

}

} // end of internal namespace
} // end of integrate namespace

#endif //INTEGRATE_LLG_ROTATION_H_
