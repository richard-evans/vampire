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
#include "random.hpp"
#include "errors.hpp"
#include "vio.hpp"

// Local cuda headers

#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"

#include "exchange_fields.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda{

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(){

#ifdef CUDA

      std::cout << "CUDA has been enabled.\n";
//      std::cout << "compiled with " << COMP << std::endl;

      // Check if there is a compatible devices
      int n_devices;
      cudaError_t error = cudaGetDeviceCount(&n_devices);

      if( error == cudaErrorNoDevice )
      {
         std::cout << "Error: CUDA is enabled but no CUDA compatible devices are available." << std::endl;
         zlog << zTs() << "Error: CUDA is enabled but no CUDA compatible devices are available." << std::endl;
         ::err::vexit();
      }
      else if ( error == cudaErrorInsufficientDriver )
      {
         std::cout     << "Error: CUDA is enabled but no CUDA drivers are incompatible. Please update drivers." << std::endl;
         zlog << zTs() << "Error: CUDA is enabled but no CUDA drivers are incompatible. Please update drivers." << std::endl;
         ::err::vexit();
      }
      else if ( error != cudaSuccess)
      {
         std::cout     << "Error: CUDA is enabled but error querying devices." << std::endl;
         zlog << zTs() << "Error: CUDA is enabled but error querying devices." << std::endl;
         ::err::vexit();
      }





      bool success = true;

      /*
       * Set the block_size according to the number of atoms
       */

      size_t _grid_size = ::atoms::num_atoms / cu::block_size + 1UL;

      if (_grid_size < cu::grid_size)
         cu::grid_size = _grid_size;

      success = success && cu::__initialize_atoms ();
      success = success && cu::__initialize_fields ();
      success = success && cu::__initialize_cells ();
      success = success && cu::__initialize_materials ();
      success = success && cu::__initialize_topology ();
      success = success && cu::__initialize_curand ();

      // Set up the exchange fields
      if( cu::exchange::initialise_exchange() != EXIT_SUCCESS)
      {
         std::cerr << "Failed to initialise exchange" << std::endl;
         success = false;
      }

      // Successful initialization
      return success;

#else
      // Default (initializtion failed)
      return false;
#endif
   }

#ifdef CUDA

   namespace internal {

      bool __initialize_atoms ()
      {
         /*
          * Allocate memory in the device and transfer the
          * spins of the atoms.
          */

         cu::atoms::x_spin_array.resize(::atoms::num_atoms);
         cu::atoms::y_spin_array.resize(::atoms::num_atoms);
         cu::atoms::z_spin_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::x_spin_array.begin(),
               ::atoms::x_spin_array.end(),
               cu::atoms::x_spin_array.begin()
               );

         thrust::copy(
               ::atoms::y_spin_array.begin(),
               ::atoms::y_spin_array.end(),
               cu::atoms::y_spin_array.begin()
               );

         thrust::copy(
               ::atoms::z_spin_array.begin(),
               ::atoms::z_spin_array.end(),
               cu::atoms::z_spin_array.begin()
               );

         /*
          * Allocate memory in the device and transfer the
          * coordinates of the atoms.
          */

         cu::atoms::x_coord_array.resize(::atoms::num_atoms);
         cu::atoms::y_coord_array.resize(::atoms::num_atoms);
         cu::atoms::z_coord_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::x_coord_array.begin(),
               ::atoms::x_coord_array.end(),
               cu::atoms::x_coord_array.begin()
               );

         thrust::copy(
               ::atoms::y_coord_array.begin(),
               ::atoms::y_coord_array.end(),
               cu::atoms::y_coord_array.begin()
               );

         thrust::copy(
               ::atoms::z_coord_array.begin(),
               ::atoms::z_coord_array.end(),
               cu::atoms::z_coord_array.begin()
               );

         /*
          * Allocate memory and send information about the types of
          * atoms
          */

         cu::atoms::type_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::type_array.begin(),
               ::atoms::type_array.end(),
               cu::atoms::type_array.begin()
               );

         /*
          * Allocate memory and pass the cell information
          */

         cu::atoms::cell_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::cell_array.begin(),
               ::atoms::cell_array.end(),
               cu::atoms::cell_array.begin()
               );

         /*
          * Allocate the memory for the unrolled spin norm array
          */

         cu::atoms::spin_norm_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::m_spin_array.begin(),
               ::atoms::m_spin_array.end(),
               cu::atoms::spin_norm_array.begin()
               );

         return true;
      }

      bool __initialize_fields ()
      {
         /*
          * Allocate memory in the device and transfer the
          * total spin field in each atom.
          */

         cu::x_total_spin_field_array.resize(::atoms::num_atoms);
         cu::y_total_spin_field_array.resize(::atoms::num_atoms);
         cu::z_total_spin_field_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::x_total_spin_field_array.begin(),
               ::atoms::x_total_spin_field_array.end(),
               cu::x_total_spin_field_array.begin()
               );

         thrust::copy(
               ::atoms::y_total_spin_field_array.begin(),
               ::atoms::y_total_spin_field_array.end(),
               cu::y_total_spin_field_array.begin()
               );

         thrust::copy(
               ::atoms::z_total_spin_field_array.begin(),
               ::atoms::z_total_spin_field_array.end(),
               cu::z_total_spin_field_array.begin()
               );

         /*
          * Allocate memory in the device and transfer the
          * total external field in each atom.
          */

         cu::x_total_external_field_array.resize(::atoms::num_atoms);
         cu::y_total_external_field_array.resize(::atoms::num_atoms);
         cu::z_total_external_field_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::x_total_external_field_array.begin(),
               ::atoms::x_total_external_field_array.end(),
               cu::x_total_external_field_array.begin()
               );

         thrust::copy(
               ::atoms::y_total_external_field_array.begin(),
               ::atoms::y_total_external_field_array.end(),
               cu::y_total_external_field_array.begin()
               );

         thrust::copy(
               ::atoms::z_total_external_field_array.begin(),
               ::atoms::z_total_external_field_array.end(),
               cu::z_total_external_field_array.begin()
               );

         /*
          * Allocate memory and transfer any existing
          * initial data for the dipolar field
          */

         cu::x_dipolar_field_array.resize(::atoms::num_atoms);
         cu::y_dipolar_field_array.resize(::atoms::num_atoms);
         cu::z_dipolar_field_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::x_dipolar_field_array.begin(),
               ::atoms::x_dipolar_field_array.end(),
               cu::x_dipolar_field_array.begin()
               );

         thrust::copy(
               ::atoms::y_dipolar_field_array.begin(),
               ::atoms::y_dipolar_field_array.end(),
               cu::y_dipolar_field_array.begin()
               );

         thrust::copy(
               ::atoms::z_dipolar_field_array.begin(),
               ::atoms::z_dipolar_field_array.end(),
               cu::z_dipolar_field_array.begin()
               );

         return true;
      }

      bool __initialize_cells ()
      {
         /*
          * Allocate memory and initialize coordinates
          */

         cu::cells::x_coord_array.resize(::cells::num_cells);
         cu::cells::y_coord_array.resize(::cells::num_cells);
         cu::cells::z_coord_array.resize(::cells::num_cells);

         thrust::copy(
               ::cells::x_coord_array.begin(),
               ::cells::x_coord_array.end(),
               cu::cells::x_coord_array.begin()
               );

         thrust::copy(
               ::cells::y_coord_array.begin(),
               ::cells::y_coord_array.end(),
               cu::cells::y_coord_array.begin()
               );

         thrust::copy(
               ::cells::z_coord_array.begin(),
               ::cells::z_coord_array.end(),
               cu::cells::z_coord_array.begin()
               );

         /*
          * Allocate memory and initialize cell magnetization
          */

         cu::cells::x_mag_array.resize(::cells::num_cells);
         cu::cells::y_mag_array.resize(::cells::num_cells);
         cu::cells::z_mag_array.resize(::cells::num_cells);

         thrust::copy(
               ::cells::x_mag_array.begin(),
               ::cells::x_mag_array.end(),
               cu::cells::x_mag_array.begin()
               );

         thrust::copy(
               ::cells::y_mag_array.begin(),
               ::cells::y_mag_array.end(),
               cu::cells::y_mag_array.begin()
               );

         thrust::copy(
               ::cells::z_mag_array.begin(),
               ::cells::z_mag_array.end(),
               cu::cells::z_mag_array.begin()
               );

         /*
          * Allocate memory and initialize cell fields
          */

         cu::cells::x_field_array.resize(::cells::num_cells);
         cu::cells::y_field_array.resize(::cells::num_cells);
         cu::cells::z_field_array.resize(::cells::num_cells);

         thrust::copy(
               ::cells::x_field_array.begin(),
               ::cells::x_field_array.end(),
               cu::cells::x_field_array.begin()
               );

         thrust::copy(
               ::cells::y_field_array.begin(),
               ::cells::y_field_array.end(),
               cu::cells::y_field_array.begin()
               );

         thrust::copy(
               ::cells::z_field_array.begin(),
               ::cells::z_field_array.end(),
               cu::cells::z_field_array.begin()
               );

         /*
          * Copy volume and number of atoms for each cell
          */

         cu::cells::volume_array.resize(::cells::num_cells);

         thrust::copy(
               ::cells::volume_array.begin(),
               ::cells::volume_array.end(),
               cu::cells::volume_array.begin()
               );

         cu::cells::num_atoms.resize(::cells::num_cells);

         thrust::copy(
               ::cells::num_atoms_in_cell.begin(),
               ::cells::num_atoms_in_cell.end(),
               cu::cells::num_atoms.begin()
               );

         return true;
      }

      bool __initialize_materials ()
      {

         /*
          * Serialize material data
          */
         size_t num_mats = ::mp::num_materials;
         thrust::host_vector<material_parameters_t> _materials(num_mats);
         for (size_t i = 0; i < num_mats; i++)
         {
            double mu_s_SI = ::mp::material[i].mu_s_SI;

            _materials[i].alpha =
               ::mp::material[i].alpha;
            _materials[i].gamma_rel =
               ::mp::material[i].gamma_rel;
            _materials[i].mu_s_si =
               mu_s_SI;
            _materials[i].i_mu_s_si =
               1.0 / mu_s_SI;
            _materials[i].k_latt =
               ::mp::material[i].Klatt_SI / mu_s_SI;
            _materials[i].sh2 =
               ::mp::material[i].sh2 / mu_s_SI;
            _materials[i].sh4 =
               ::mp::material[i].sh4 / mu_s_SI;
            _materials[i].sh6 =
               ::mp::material[i].sh6 / mu_s_SI;
            _materials[i].ku =
               ::mp::material[i].Ku;
            _materials[i].anisotropy_unit_x =
               ::mp::material[i].UniaxialAnisotropyUnitVector[0];
            _materials[i].anisotropy_unit_y =
               ::mp::material[i].UniaxialAnisotropyUnitVector[1];
            _materials[i].anisotropy_unit_z =
               ::mp::material[i].UniaxialAnisotropyUnitVector[2];
            _materials[i].applied_field_strength =
               ::mp::material[i].applied_field_strength;
            _materials[i].applied_field_unit_x =
               ::mp::material[i].applied_field_unit_vector[0];
            _materials[i].applied_field_unit_y =
               ::mp::material[i].applied_field_unit_vector[1];
            _materials[i].applied_field_unit_z =
               ::mp::material[i].applied_field_unit_vector[2];
            _materials[i].Kc1_SI =
               ::mp::material[i].Kc1_SI;
            _materials[i].temperature =
               ::mp::material[i].temperature;
            _materials[i].temperature_rescaling_alpha =
               ::mp::material[i].temperature_rescaling_alpha;
            _materials[i].temperature_rescaling_Tc =
               ::mp::material[i].temperature_rescaling_Tc;
            _materials[i].H_th_sigma =
               ::mp::material[i].H_th_sigma;
         }

         /*
          * Allocate memory and send information about the materials
          */
         cu::mp::materials.resize(num_mats);
         thrust::copy(
            _materials.begin(),
            _materials.end(),
            cu::mp::materials.begin()
            );

         return true;
      }

      bool __initialize_topology ()
      {
         /*
          * Send the information for limits and neighbors up to the
          * device.
          *

         // Resize and set all values to 0
         cu::atoms::limits.assign(::atoms::num_atoms + 1UL, 0);
         cu::atoms::neighbours.resize(::atoms::total_num_neighbours);

         thrust::copy(
               ::atoms::neighbour_list_end_index.begin(),
               ::atoms::neighbour_list_end_index.end(),
               cu::atoms::limits.begin() + 1UL
               );

         *
          * Transform the limits to be one pased the last element
          * in the neighbors list.
          *
         thrust::transform(
               cu::atoms::limits.begin(),
               cu::atoms::limits.end(),
               cu::atoms::limits.begin(),
               cu::plusone_functor()
               );

         thrust::copy(
               ::atoms::neighbour_list_array.begin(),
               ::atoms::neighbour_list_array.end(),
               cu::atoms::neighbours.begin()
               );
         */

         // Transfer the row ptrs and col indices to the device
         std::vector<int> limits_h( ::atoms::num_atoms + 1, 0);
         for( int atom = 0; atom < ::atoms::num_atoms; atom++)
            limits_h[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;

         cu::atoms::limits.resize( ::atoms::num_atoms + 1);
         cu::atoms::neighbours.resize( ::atoms::neighbour_list_array.size() );


         thrust::copy(
               limits_h.begin(),
               limits_h.end(),
               cu::atoms::limits.begin()
               );

         thrust::copy(
               ::atoms::neighbour_list_array.begin(),
               ::atoms::neighbour_list_array.end(),
               cu::atoms::neighbours.begin()
               );



         return true;
      }

      bool __initialize_curand ()
      {
         cudaMalloc (
               (void **) &cu::d_rand_state,
               cu::grid_size * cu::block_size * sizeof(curandState));

         check_cuda_errors (__FILE__, __LINE__);

         cu::init_rng <<< cu::grid_size, cu::block_size >>> (
               cu::d_rand_state, ::mtrandom::integration_seed);

         check_cuda_errors (__FILE__, __LINE__);

         return true;
      }
   }

#endif

} // end of namespace vcuda
