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

      bool success = true;

      success = success || __initialize_atoms ();
      success = success || __initialize_fields ();
      success = success || __initialize_cells ();
      success = success || __initialize_materials ();
      success = success || __initialize_topology ();

      // send topology information

      // Successful initialization
      return success;

#endif

      // Default (initializtion failed)
      return false;

   }

   bool __initialize_atoms ()
   {
      /*
       * Allocate memory in the device and transfer the
       * spins of the atoms.
       */

      cu::x_spin_array.resize(atoms::num_atoms);
      cu::y_spin_array.resize(atoms::num_atoms);
      cu::z_spin_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_spin_array.begin(),
            atoms::x_spin_array.end(),
            cu::x_spin_array.begin()
            );

      thrust::copy(
            atoms::y_spin_array.begin(),
            atoms::y_spin_array.end(),
            cu::y_spin_array.begin()
            );

      thrust::copy(
            atoms::z_spin_array.begin(),
            atoms::z_spin_array.end(),
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
            atoms::x_coord_array.end(),
            cu::x_coord_array.begin()
            );

      thrust::copy(
            atoms::y_coord_array.begin(),
            atoms::y_coord_array.end(),
            cu::y_coord_array.begin()
            );

      thrust::copy(
            atoms::z_coord_array.begin(),
            atoms::z_coord_array.end(),
            cu::z_coord_array.begin()
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
       * Allocate memory and pass the cell information
       */

      cu::cell_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::cell_array.begin(),
            atoms::cell_array.end(),
            cu::cell_array.begin()
            );
   }

   bool __initialize_fields ()
   {
      /*
       * Allocate memory in the device and transfer the
       * total spin field in each atom.
       */

      cu::x_total_spin_field_array.resize(atoms::num_atoms);
      cu::y_total_spin_field_array.resize(atoms::num_atoms);
      cu::z_total_spin_field_array.resize(atoms::num_atoms);

      thrust::copy(
            atoms::x_total_spin_field_array.begin(),
            atoms::x_total_spin_field_array.end(),
            cu::x_total_spin_field_array.begin()
            );

      thrust::copy(
            atoms::y_total_spin_field_array.begin(),
            atoms::y_total_spin_field_array.end(),
            cu::y_total_spin_field_array.begin()
            );

      thrust::copy(
            atoms::z_total_spin_field_array.begin(),
            atoms::z_total_spin_field_array.end(),
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
            atoms::x_total_external_field_array.end(),
            cu::x_total_external_field_array.begin()
            );

      thrust::copy(
            atoms::y_total_external_field_array.begin(),
            atoms::y_total_external_field_array.end(),
            cu::y_total_external_field_array.begin()
            );

      thrust::copy(
            atoms::z_total_external_field_array.begin(),
            atoms::z_total_external_field_array.end(),
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
            atoms::x_dipolar_field_array.end(),
            cu::x_dipolar_field_array.begin()
            );

      thrust::copy(
            atoms::y_dipolar_field_array.begin(),
            atoms::y_dipolar_field_array.end(),
            cu::y_dipolar_field_array.begin()
            );

      thrust::copy(
            atoms::z_dipolar_field_array.begin(),
            atoms::z_dipolar_field_array.end(),
            cu::z_dipolar_field_array.begin()
            );
   }

   bool __initialize_cells ()
   {
      /*
       * Allocate memory and initialize coordinates
       */

      cu::cell_x_coord_array.resize(cells::num_cells);
      cu::cell_y_coord_array.resize(cells::num_cells);
      cu::cell_z_coord_array.resize(cells::num_cells);

      thrust::copy(
            cells::x_coord_array.begin(),
            cells::x_coord_array.end(),
            cu::cell_x_coord_array.begin()
            );

      thrust::copy(
            cells::y_coord_array.begin(),
            cells::y_coord_array.end(),
            cu::cell_y_coord_array.begin()
            );

      thrust::copy(
            cells::z_coord_array.begin(),
            cells::z_coord_array.end(),
            cu::cell_z_coord_array.begin()
            );

      /*
       * Allocate memory and initialize cell magnetization
       */

      cu::cell_x_mag_array.resize(cells::num_cells);
      cu::cell_y_mag_array.resize(cells::num_cells);
      cu::cell_z_mag_array.resize(cells::num_cells);

      thrust::copy(
            cells::x_mag_array.begin(),
            cells::x_mag_array.end(),
            cu::cell_x_mag_array.begin()
            );

      thrust::copy(
            cells::y_mag_array.begin(),
            cells::y_mag_array.end(),
            cu::cell_y_mag_array.begin()
            );

      thrust::copy(
            cells::z_mag_array.begin(),
            cells::z_mag_array.end(),
            cu::cell_z_mag_array.begin()
            );

      /*
       * Copy volume and number of atoms for each cell
       */

      cu::cell_volume_array.resize(cells::num_cells);

      thrust::copy(
            cells::volume_array.begin(),
            cells::volume_array.end(),
            cu::cell_volume_array.begin()
            );

      cu::cell_num_atoms.resize(cells::num_cells);

      thrust::copy(
            cells::num_atoms_in_cell.begin(),
            cells::num_atoms_in_cell.end(),
            cu::cell_num_atoms.begin()
            );
   }

   bool __initialize_materials ()
   {
      /*
       * Allocate memory and send information about the materias
       */

      cu::materials.resize(mp::num_materials);

      thrust::copy(
            mp::material.begin(),
            mp::material.end(),
            cu::materials.begin()
            );

   }

   bool __initialize_topology ()
   {
      /*
       * Send the information for limits and neighbors up to the
       * device.
       */

      cu::limits.resize(atoms::num_atoms);
      cu::neighbours.resize(atoms::total_num_neighbours);

      thrust::copy(
            atoms::neighbour_list_end_index.begin(),
            atoms::neighbour_list_end_index.end(),
            cu::limits.begin()
            );

      /*
       * Transform the limits to be one pased the last element
       * in the neighbors list.
       */
      thrust::transform(
            cu::limits.begin(),
            cu::limits.end(),
            cu::limits.begin(),
            cu::plusone_functor()
            );

      thrust::copy(
            atoms::neighbour_list_array.begin(),
            atoms::neighbour_list_array.end(),
            cu::neighbours.begin()
            );

   }

} // end of namespace cuda
