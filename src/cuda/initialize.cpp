//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

#include <vector>

// Thrust library headers

#include <thrust/device_vector.h>
#include <thrust/copy.h>

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "internal.hpp"

namespace cu = cuda::internal;

namespace cuda{

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

#ifdef CUDA

      // send the spin information

      /*
       * Allocate memory in the device and transfer the
       * spins of the atoms.
       */

      cu::x_spin_array.resize(atoms::num_atoms);
      cu::y_spin_array.resize(atoms::num_atoms);
      cu::z_spin_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_spin_array.begin(),
            atoms::x_spin_array.end()
            cu::x_spin_array.begin()
            );

      thrust::copy(
            atoms::y_spin_array.begin(),
            atoms::y_spin_array.end()
            cu::y_spin_array.begin()
            );

      thrust::copy(
            cu::atoms::z_spin_array.begin(),
            cu::atoms::z_spin_array.end()
            cu::z_spin_array.begin()
            );

      /*
       * Allocate memory in the device and transfer the
       * coordinates of the atoms.
       */

      cu::x_coord_array.resize(atoms::num_atoms);
      cu::y_coord_array.resize(atoms::num_atoms);
      cu::z_coord_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_coord_array.begin(),
            atoms::x_coord_array.end()
            cu::x_coord_array.begin()
            );

      thrust::copy(
            atoms::y_coord_array.begin(),
            atoms::y_coord_array.end()
            cu::y_coord_array.begin()
            );

      thrust::copy(
            atoms::z_coord_array.begin(),
            atoms::z_coord_array.end()
            cu::z_coord_array.begin()
            );

      /*
       * Allocate memory in the device and transfer the
       * total spin field in each atom.
       */

      cu::x_total_spin_field_array.resize(atoms::num_atoms);
      cu::y_total_spin_field_array.resize(atoms::num_atoms);
      cu::z_total_spin_field_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_total_spin_field_array.begin(),
            atoms::x_total_spin_field_array.end()
            cu::x_total_spin_field_array.begin()
            );

      thrust::copy(
            atoms::y_total_spin_field_array.begin(),
            atoms::y_total_spin_field_array.end()
            cu::y_total_spin_field_array.begin()
            );

      thrust::copy(
            atoms::z_total_spin_field_array.begin(),
            atoms::z_total_spin_field_array.end()
            cu::z_total_spin_field_array.begin()
            );

      /*
       * Allocate memory in the device and transfer the
       * total external field in each atom.
       */

      cu::x_total_external_field_array.resize(atoms::num_atoms);
      cu::y_total_external_field_array.resize(atoms::num_atoms);
      cu::z_total_external_field_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_total_external_field_array.begin(),
            atoms::x_total_external_field_array.end()
            cu::x_total_external_field_array.begin()
            );

      thrust::copy(
            atoms::y_total_external_field_array.begin(),
            atoms::y_total_external_field_array.end()
            cu::y_total_external_field_array.begin()
            );

      thrust::copy(
            atoms::z_total_external_field_array.begin(),
            atoms::z_total_external_field_array.end()
            cu::z_total_external_field_array.begin()
            );

      /*
       * Allocate memory and transfer any existing
       * initial data for the dipolar field
       */

      cu::x_dipolar_field_array.resize(atoms::num_atoms);
      cu::y_dipolar_field_array.resize(atoms::num_atoms);
      cu::z_dipolar_field_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_dipolar_field_array.begin(),
            atoms::x_dipolar_field_array.end()
            cu::x_dipolar_field_array.begin()
            );

      thrust::copy(
            atoms::y_dipolar_field_array.begin(),
            atoms::y_dipolar_field_array.end()
            cu::y_dipolar_field_array.begin()
            );

      thrust::copy(
            atoms::z_dipolar_field_array.begin(),
            atoms::z_dipolar_field_array.end()
            cu::z_dipolar_field_array.begin()
            );

      /*
       * Allocate memory and send information about the types of
       * atoms
       */

      cu::type_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::type_array.begin(),
            atoms::type_array.end(),
            cu::type_array.begin()
            );

      /*
       * Allocate memory and send information about the materias
       */

      cu::materias.resize(mp::num_materials);

      thrust::copy(
            mp::materials.begin(),
            mp::materials.end(),
            cu::materials.begin()
            );

      // send the material information
      // send topology information
      // send the macro-cell information

      // Successful initialization
      return true;

#endif

      // Default (initializtion failed)
      return false;

   }

} // end of namespace cuda
