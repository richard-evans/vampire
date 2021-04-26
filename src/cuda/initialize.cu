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
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "cuda.hpp"
#include "errors.hpp"
#include "dipole.hpp"
#include "gpu.hpp"
#include "random.hpp"
#include "stats.hpp"
#include "typedefs.hpp"
#include "vio.hpp"

// Local cuda headers

#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"

#include "exchange_fields.hpp"
#include "statistics.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda{

   //-------------------------------------------------------------------------------
   // Function to initialize GPU data
   //-------------------------------------------------------------------------------
   bool initialize(bool cpu_stats){

#ifdef CUDA

      std::cout << "CUDA has been enabled in ";
      zlog << zTs() << "CUDA has been enabled in ";
      #ifdef CUDA_DP
         std::cout << "double precision mode" << std::endl;
         zlog << "double precision mode" << std::endl;
      #else
         std::cout << "single precision mode" << std::endl;
         zlog << "single precision mode" << std::endl;
      #endif

      // set internal cpu statistics flag
      vcuda::internal::stats::use_cpu = cpu_stats;

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

      // Set cuda device if specified by user
      if(gpu::device != -1 && gpu::device < n_devices){
         zlog << zTs() << "Setting CUDA device to " << gpu::device << std::endl;
         cudaError_t error = cudaSetDevice(gpu::device);
         if( error == cudaErrorSetOnActiveProcess )
         {
            std::cerr     << "Error: CUDA is unable to set active process to device " << gpu::device << std::endl;
            zlog << zTs() << "Error: CUDA is unable to set active process to device " << gpu::device << std::endl;
            ::err::vexit();
         }
         else if ( error == cudaErrorInvalidDevice )
         {
            std::cerr     << "Error: CUDA is requesting device " << gpu::device << " which is an invalid device." << std::endl;
            zlog << zTs() << "Error: CUDA is requesting device " << gpu::device << " which is an invalid device." << std::endl;
            ::err::vexit();
         }
         else if ( error != cudaSuccess)
         {
            std::cerr     << "Error: CUDA unable to set device to " << gpu::device << std::endl;
            zlog << zTs() << "Error: CUDA unable to set device to " << gpu::device << std::endl;
            ::err::vexit();
         }
      }

      bool success = true;

      /*
       * Set the block_size according to the number of atoms
       */

      size_t _grid_size = ( (::atoms::num_atoms/2) / cu::block_size) + 1;

      //std::cerr << "Natoms = " << ::atoms::num_atoms << "\n";
      //std::cerr << "Block size = " << cu::block_size << "\n";
      //std::cerr << "grid1 = " << _grid_size << "\n";
      //std::cerr << "grid2 = " << cu::grid_size << "\n";

      // I dont think this does anything
      //if (_grid_size < cu::grid_size)
      //   cu::grid_size = _grid_size;

      cu::grid_size = _grid_size;
      //std::cerr << "grid2 = " << cu::grid_size << std::endl;

      zlog << zTs() << "Using cuda version with block size " << cu::block_size << " and grid size " << cu::grid_size << std::endl;

      success = success && cu::__initialize_atoms ();
      success = success && cu::__initialize_fields ();
      success = success && cu::__initialize_cells ();
      success = success && cu::__initialize_materials ();
      success = success && cu::__initialize_topology ();
      success = success && cu::__initialize_curand ();
      success = success && cu::__initialize_stats ();

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
   

   bool initialize_dipole(){
#ifdef CUDA

      bool success = true;

      // Initialise dipole
      if( cu::__initialize_dipole() != true)
      {
         std::cerr << "Failed to initialise dipole" << std::endl;
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
         /*
         cu::atoms::x_spin_array.resize(::atoms::num_atoms);
         cu::atoms::y_spin_array.resize(::atoms::num_atoms);
         cu::atoms::z_spin_array.resize(::atoms::num_atoms);
         */
         cudaMalloc((void**)&cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu_real_t));
         /* Need to be careful here
         The device code can use SP or SP,
         but the host code seems to rely exclusively on DP */
         cudaMemcpy(cu::atoms::d_x_spin, ::atoms::x_spin_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::atoms::d_y_spin, ::atoms::y_spin_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::atoms::d_z_spin, ::atoms::z_spin_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*thrust::copy(
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
         */
         /*
          * Allocate memory in the device and transfer the
          * coordinates of the atoms.
          */
         /*
         cu::atoms::x_coord_array.resize(::atoms::num_atoms);
         cu::atoms::y_coord_array.resize(::atoms::num_atoms);
         cu::atoms::z_coord_array.resize(::atoms::num_atoms);
         */

         cudaMalloc((void**)&cu::atoms::d_x_coord, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::atoms::d_y_coord, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::atoms::d_z_coord, ::atoms::num_atoms * sizeof(cu_real_t));

         cudaMemcpy(cu::atoms::d_x_coord, ::atoms::x_coord_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::atoms::d_y_coord, ::atoms::y_coord_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::atoms::d_z_coord, ::atoms::z_coord_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);


         /*thrust::copy(
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
         */
         /*
          * Allocate memory and send information about the types of
          * atoms
          */

         //cu::atoms::type_array.resize(::atoms::num_atoms);

         cudaMalloc((void**)&cu::atoms::d_materials, ::atoms::num_atoms * sizeof(int));
         cudaMemcpy(cu::atoms::d_materials, ::atoms::type_array.data(), ::atoms::num_atoms * sizeof(int), cudaMemcpyHostToDevice);

         /*thrust::copy(
               ::atoms::type_array.begin(),
               ::atoms::type_array.end(),
               cu::atoms::type_array.begin()
               );
         */
         /*
          * Allocate memory and pass the cell information
          */

         cudaMalloc((void**)&cu::atoms::d_cells, ::atoms::num_atoms * sizeof(int));
         cudaMemcpy(cu::atoms::d_cells, ::atoms::cell_array.data(), ::atoms::num_atoms * sizeof(int), cudaMemcpyHostToDevice);


         /*cu::atoms::cell_array.resize(::atoms::num_atoms);

         thrust::copy(
               ::atoms::cell_array.begin(),
               ::atoms::cell_array.end(),
               cu::atoms::cell_array.begin()
               );
         */
         /*
          * Allocate the memory for the unrolled spin norm array
          */

         // This is actually used in Thrust algorithms in statistics.cu
         // Leave it for now
         //cu::atoms::spin_norm_array.resize(::atoms::num_atoms);

         //thrust::copy(
         //      ::atoms::m_spin_array.begin(),
         //      ::atoms::m_spin_array.end(),
         //      cu::atoms::spin_norm_array.begin()
         //      );

         return true;
      }

      bool __initialize_fields ()
      {
         /*
          * Allocate memory in the device and transfer the
          * total spin field in each atom.
          */

         cudaMalloc((void**)&cu::d_x_spin_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_y_spin_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_z_spin_field, ::atoms::num_atoms * sizeof(cu_real_t));

         cudaMemcpy(cu::d_x_spin_field, ::atoms::x_total_spin_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_y_spin_field, ::atoms::y_total_spin_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_z_spin_field, ::atoms::z_total_spin_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
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
         */

         /*
          * Allocate memory in the device and transfer the
          * total external field in each atom.
          */

         cudaMalloc((void**)&cu::d_x_external_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_y_external_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_z_external_field, ::atoms::num_atoms * sizeof(cu_real_t));

         cudaMemcpy(cu::d_x_external_field, ::atoms::x_total_external_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_y_external_field, ::atoms::y_total_external_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_z_external_field, ::atoms::z_total_external_field_array.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);


         /*cu::x_total_external_field_array.resize(::atoms::num_atoms);
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
         */
         /*
          * Allocate memory and transfer any existing
          * initial data for the dipolar field
          */

         cudaMalloc((void**)&cu::d_x_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_y_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_z_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));

         cudaMemcpy(cu::d_x_dip_field, ::dipole::atom_dipolar_field_array_x.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_y_dip_field, ::dipole::atom_dipolar_field_array_y.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_z_dip_field, ::dipole::atom_dipolar_field_array_z.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         cudaMalloc((void**)&cu::d_x_mu0H_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_y_mu0H_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::d_z_mu0H_dip_field, ::atoms::num_atoms * sizeof(cu_real_t));

         cudaMemcpy(cu::d_x_mu0H_dip_field, ::dipole::atom_mu0demag_field_array_x.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_y_mu0H_dip_field, ::dipole::atom_mu0demag_field_array_x.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::d_z_mu0H_dip_field, ::dipole::atom_mu0demag_field_array_x.data(), ::atoms::num_atoms * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
         cu::x_dipolar_field_array.resize(::atoms::num_atoms);
         cu::y_dipolar_field_array.resize(::atoms::num_atoms);
         cu::z_dipolar_field_array.resize(::atoms::num_atoms);

         thrust::copy(::dipole::atom_dipolar_field_array_x.begin(),::dipole::atom_dipolar_field_array_x.end(), cu::x_dipolar_field_array.begin());
         thrust::copy(::dipole::atom_dipolar_field_array_y.begin(),::dipole::atom_dipolar_field_array_y.end(), cu::y_dipolar_field_array.begin());
         thrust::copy(::dipole::atom_dipolar_field_array_z.begin(),::dipole::atom_dipolar_field_array_z.end(), cu::z_dipolar_field_array.begin());
         */
         return true;
      }

      bool __initialize_cells ()
      {
         /*
          * Allocate memory and initialize coordinates
          */

         cudaMalloc((void**)&cu::cells::d_x_coord, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_y_coord, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_z_coord, ::cells::num_cells * sizeof(cu_real_t));


         /*cu::cells::x_coord_array.resize(::cells::num_cells);
         cu::cells::y_coord_array.resize(::cells::num_cells);
         cu::cells::z_coord_array.resize(::cells::num_cells);
         */
         // unroll 4N array to N
         std::vector<double> pos(::cells::num_cells,0.0);
         for(int cell = 0; cell < pos.size(); cell++) pos[cell] = ::cells::pos_and_mom_array[4*cell+0]; // x

         //thrust::copy(pos.begin(), pos.end(), cu::cells::x_coord_array.begin());
         cudaMemcpy(cu::cells::d_x_coord, pos.data(), ::cells::num_cells, cudaMemcpyHostToDevice);

         // unroll 4N array to N
         for(int cell = 0; cell < pos.size(); cell++) pos[cell] = ::cells::pos_and_mom_array[4*cell+1]; // y

         //thrust::copy(pos.begin(), pos.end(), cu::cells::y_coord_array.begin());
         cudaMemcpy(cu::cells::d_y_coord, pos.data(), ::cells::num_cells, cudaMemcpyHostToDevice);

         // unroll 4N array to N
         for(int cell = 0; cell < pos.size(); cell++) pos[cell] = ::cells::pos_and_mom_array[4*cell+2]; // z

         //thrust::copy(pos.begin(), pos.end(), cu::cells::z_coord_array.begin());
         cudaMemcpy(cu::cells::d_z_coord, pos.data(), ::cells::num_cells, cudaMemcpyHostToDevice);

         //-----------------------------------------------------
         // Allocate memory and initialize cell magnetization
         //-----------------------------------------------------

         cudaMalloc((void**)&cu::cells::d_x_mag, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_y_mag, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_z_mag, ::cells::num_cells * sizeof(cu_real_t));

         cudaMemcpy(cu::cells::d_x_mag, ::cells::mag_array_x.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_y_mag, ::cells::mag_array_y.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_z_mag, ::cells::mag_array_z.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
         cu::cells::x_mag_array.resize(::cells::num_cells);
         cu::cells::y_mag_array.resize(::cells::num_cells);
         cu::cells::z_mag_array.resize(::cells::num_cells);

         thrust::copy(::cells::mag_array_x.begin(), ::cells::mag_array_x.end(), cu::cells::x_mag_array.begin());
         thrust::copy(::cells::mag_array_y.begin(), ::cells::mag_array_y.end(), cu::cells::y_mag_array.begin());
         thrust::copy(::cells::mag_array_z.begin(), ::cells::mag_array_z.end(), cu::cells::z_mag_array.begin());
         */

         //----------------------------------------------
         // Allocate memory and initialize cell fields
         //----------------------------------------------

         cudaMalloc((void**)&cu::cells::d_x_cell_field, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_y_cell_field, ::cells::num_cells * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_z_cell_field, ::cells::num_cells * sizeof(cu_real_t));

         cudaMemcpy(cu::cells::d_x_cell_field, ::cells::field_array_x.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_y_cell_field, ::cells::field_array_x.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_z_cell_field, ::cells::field_array_x.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
         cu::cells::x_field_array.resize(::cells::num_cells);
         cu::cells::y_field_array.resize(::cells::num_cells);
         cu::cells::z_field_array.resize(::cells::num_cells);
         thrust::copy(::cells::field_array_x.begin(), ::cells::field_array_x.end(), cu::cells::x_field_array.begin());
         thrust::copy(::cells::field_array_y.begin(), ::cells::field_array_y.end(), cu::cells::y_field_array.begin());
         thrust::copy(::cells::field_array_z.begin(), ::cells::field_array_z.end(), cu::cells::z_field_array.begin());
         */

         /*
          * Copy volume and number of atoms for each cell
          */

         cudaMalloc((void**)&cu::cells::d_volume, ::cells::num_cells * sizeof(cu_real_t));
         cudaMemcpy(cu::cells::d_volume, ::cells::volume_array.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
         cu::cells::volume_array.resize(::cells::num_cells);

         thrust::copy(
               ::cells::volume_array.begin(),
               ::cells::volume_array.end(),
               cu::cells::volume_array.begin()
               );
         */

         cudaMalloc((void**)&cu::cells::d_num_atoms, ::cells::num_cells * sizeof(cu_real_t));
         cudaMemcpy(cu::cells::d_num_atoms, ::cells::num_atoms_in_cell.data(), ::cells::num_cells * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         /*
         cu::cells::num_atoms.resize(::cells::num_cells);

         thrust::copy(
               ::cells::num_atoms_in_cell.begin(),
               ::cells::num_atoms_in_cell.end(),
               cu::cells::num_atoms.begin()
               );
         */
         return true;
      }
      
      bool __initialize_dipole(){
         
         std::vector<double> tensor_xx = ::dipole::get_tensor_1D_xx();
         std::vector<double> tensor_xy = ::dipole::get_tensor_1D_xy();
         std::vector<double> tensor_xz = ::dipole::get_tensor_1D_xz();
         std::vector<double> tensor_yy = ::dipole::get_tensor_1D_yy();
         std::vector<double> tensor_yz = ::dipole::get_tensor_1D_yz();
         std::vector<double> tensor_zz = ::dipole::get_tensor_1D_zz();

         cudaMalloc((void**)&cu::cells::d_tensor_xx, tensor_xx.size() * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_tensor_xy, tensor_xy.size() * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_tensor_xz, tensor_xz.size() * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_tensor_yy, tensor_yy.size() * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_tensor_yz, tensor_yz.size() * sizeof(cu_real_t));
         cudaMalloc((void**)&cu::cells::d_tensor_zz, tensor_zz.size() * sizeof(cu_real_t));
         
         check_cuda_errors(__FILE__,__LINE__);
         

         cudaMemcpy(cu::cells::d_tensor_xx, tensor_xx.data(), tensor_xx.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_tensor_xy, tensor_xy.data(), tensor_xy.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_tensor_xz, tensor_xz.data(), tensor_xz.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_tensor_yy, tensor_yy.data(), tensor_yy.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_tensor_yz, tensor_yz.data(), tensor_yz.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::cells::d_tensor_zz, tensor_zz.data(), tensor_zz.size() * sizeof(cu_real_t), cudaMemcpyHostToDevice);

         check_device_memory(__FILE__,__LINE__);
         
         return true;
      }

      bool __initialize_materials ()
      {

         /*
          * Serialize material data
          */
         size_t num_mats = ::mp::num_materials;
         //thrust::host_vector<material_parameters_t> _materials(num_mats);
         std::vector<material_parameters_t> _materials(num_mats);
         for (size_t i = 0; i < num_mats; i++)
         {
            double mu_s_SI = ::mp::material[i].mu_s_SI;

            double ku2 = ::anisotropy::get_ku2(i); // second order uniaxial anisotropy constant (Ku1)
            double ku4 = ::anisotropy::get_ku4(i); // fourth order uniaxial anisotropy constant (Ku2)
            double ku6 = ::anisotropy::get_ku6(i); // sixth order uniaxial anisotropy constant  (Ku3)
            double kc4 = ::anisotropy::get_kc4(i); // fourth order cubic anisotropy constant (Kc1)
            //double kc6 = ::anisotropy::get_kc6(i); // sixth order cubic anisotropy constant (Kc2)

            std::vector<double> ku_vector = ::anisotropy::get_ku_vector(i); // unit vector defining axis for uniaxial anisotropy

            _materials[i].alpha     = ::mp::material[i].alpha;
            _materials[i].gamma_rel = ::mp::material[i].gamma_rel;
            _materials[i].mu_s_si   = mu_s_SI;
            _materials[i].i_mu_s_si = 1.0 / mu_s_SI;
            _materials[i].k_latt = 0.0; //::mp::material[i].Klatt_SI / mu_s_SI;
            _materials[i].sh2 = ku2;// / mu_s_SI;
            _materials[i].sh4 = ku4;
            _materials[i].sh6 = ku6;
            _materials[i].anisotropy_unit_x = ku_vector[0];
            _materials[i].anisotropy_unit_y = ku_vector[1];
            _materials[i].anisotropy_unit_z = ku_vector[2];
            _materials[i].applied_field_strength = ::mp::material[i].applied_field_strength;
            _materials[i].applied_field_unit_x   = ::mp::material[i].applied_field_unit_vector[0];
            _materials[i].applied_field_unit_y   = ::mp::material[i].applied_field_unit_vector[1];
            _materials[i].applied_field_unit_z   = ::mp::material[i].applied_field_unit_vector[2];
            _materials[i].kc4 = kc4;
            _materials[i].temperature = ::mp::material[i].temperature;
            _materials[i].temperature_rescaling_alpha = ::mp::material[i].temperature_rescaling_alpha;
            _materials[i].temperature_rescaling_Tc    = ::mp::material[i].temperature_rescaling_Tc;
            _materials[i].H_th_sigma = ::mp::material[i].H_th_sigma;

         }

         //std::vector<double> ks_tensor = ::anisotropy::get_ku_vector(i);

         /*
          * Allocate memory and send information about the materials
          */
          cudaMalloc((void**)&cu::mp::d_material_params, num_mats * sizeof(material_parameters_t));
          cudaMemcpy(cu::mp::d_material_params, _materials.data(), num_mats * sizeof(material_parameters_t), cudaMemcpyHostToDevice);

         /*
         cu::mp::materials.resize(num_mats);
         thrust::copy(
            _materials.begin(),
            _materials.end(),
            cu::mp::materials.begin()
            );
         */
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

         cudaMalloc((void**)&cu::atoms::d_limits, (::atoms::num_atoms + 1) * sizeof(int));
         cudaMalloc((void**)&cu::atoms::d_neighbours, ::atoms::neighbour_list_array.size() * sizeof(int));

         cudaMemcpy(cu::atoms::d_limits, limits_h.data(), (::atoms::num_atoms + 1) * sizeof(int), cudaMemcpyHostToDevice);
         cudaMemcpy(cu::atoms::d_neighbours, ::atoms::neighbour_list_array.data(), ::atoms::neighbour_list_array.size() * sizeof(int), cudaMemcpyHostToDevice);

         /*
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

         */


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

      bool __initialize_stats ()
      {

         return true;

      }

      /**
       * Inits the random number generator states in the device, one per thread
       */
      __global__ void init_rng (curandState * states, int seed)
      {
         int tid = blockIdx.x * blockDim.x + threadIdx.x;
         curand_init (seed, tid, 0, &states[tid]);
      }
   }

#endif

} // end of namespace vcuda
