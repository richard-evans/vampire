//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "exchange_fields.hpp"
#include "internal.hpp"
#include "llg_heun.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda{

   //--------------------------------------------------------------------------
   // Function to perform a single heun step
   //--------------------------------------------------------------------------
   void llg_heun(){

#ifdef CUDA

      if (!cu::llg::initialized) cu::llg::__llg_init();
      cu::llg::__llg_step();

#endif

      return;
   }

#ifdef CUDA

   namespace internal {

      namespace llg
      {
         bool initialized(false);
         RealArray x_spin_buffer_array(0UL);
         RealArray y_spin_buffer_array(0UL);
         RealArray z_spin_buffer_array(0UL);
         HeunParametersArray heun_parameters(0UL);

         void __llg_init ()
         {
            /*
             * Reserve space for the buffers
             */
            cu::llg::x_spin_buffer_array.resize(::atoms::num_atoms);
            cu::llg::y_spin_buffer_array.resize(::atoms::num_atoms);
            cu::llg::z_spin_buffer_array.resize(::atoms::num_atoms);

            /*
             * Initialize heun parameters
             */
            size_t num_mats = ::mp::num_materials;
            thrust::host_vector<heun_parameters_t> _parameters(num_mats);

            for (size_t i = 0; i < num_mats; i++)
            {
               double alpha = ::mp::material[i].alpha;
               double gamma = ::mp::material[i].gamma_rel;
               _parameters[i].prefactor = gamma / (1.0 + alpha * alpha);
               /*
                * lambda is alpha (LOL)
                */
               _parameters[i].lambda_times_prefactor =
                  gamma * alpha / (1.0 + alpha * alpha);

#ifdef CUDA_DEBUG
               std::cout << "Heun parameters: "
                   << _parameters[i].prefactor << " "
                   << _parameters[i].lambda_times_prefactor << std::endl;
#endif

            }

            cu::llg::heun_parameters.resize(num_mats);
            thrust::copy(
                  _parameters.begin(),
                  _parameters.end(),
                  cu::llg::heun_parameters.begin()
                  );
         }

         void __llg_step ()
         {
            thrust::copy (
                  cu::atoms::x_spin_array.begin(),
                  cu::atoms::x_spin_array.end(),
                  cu::llg::x_spin_buffer_array.begin()
                  );
            thrust::copy (
                  cu::atoms::y_spin_array.begin(),
                  cu::atoms::y_spin_array.end(),
                  cu::llg::y_spin_buffer_array.begin()
                  );
            thrust::copy (
                  cu::atoms::z_spin_array.begin(),
                  cu::atoms::z_spin_array.end(),
                  cu::llg::z_spin_buffer_array.begin()
                  );

#ifdef CUDA_DEBUG
            std::cout << cu::atoms::x_spin_array[0] << " "
                      << cu::atoms::y_spin_array[0] << " "
                      << cu::atoms::z_spin_array[0] << std::endl;
#endif

            double * d_x_spin = thrust::raw_pointer_cast(
                  cu::atoms::x_spin_array.data());
            double * d_y_spin = thrust::raw_pointer_cast(
                  cu::atoms::y_spin_array.data());
            double * d_z_spin = thrust::raw_pointer_cast(
                  cu::atoms::z_spin_array.data());

            int * d_materials =
               thrust::raw_pointer_cast(cu::atoms::type_array.data());

            cu::heun_parameters_t * d_heun_params =
               thrust::raw_pointer_cast (cu::llg::heun_parameters.data());

            double * d_x_spin_field = thrust::raw_pointer_cast(
                  cu::x_total_spin_field_array.data());
            double * d_y_spin_field = thrust::raw_pointer_cast(
                  cu::y_total_spin_field_array.data());
            double * d_z_spin_field = thrust::raw_pointer_cast(
                  cu::z_total_spin_field_array.data());

            double * d_x_external_field = thrust::raw_pointer_cast(
                  cu::x_total_external_field_array.data());
            double * d_y_external_field = thrust::raw_pointer_cast(
                  cu::y_total_external_field_array.data());
            double * d_z_external_field = thrust::raw_pointer_cast(
                  cu::z_total_external_field_array.data());

            double * d_x_spin_buffer = thrust::raw_pointer_cast(
                  cu::llg::x_spin_buffer_array.data());
            double * d_y_spin_buffer = thrust::raw_pointer_cast(
                  cu::llg::y_spin_buffer_array.data());
            double * d_z_spin_buffer = thrust::raw_pointer_cast(
                  cu::llg::z_spin_buffer_array.data());

            cu::update_spin_fields ();
            cu::update_external_fields ();

            check_cuda_errors (__FILE__, __LINE__);

            cu::llg::llg_heun_step <<< cu::grid_size, cu::block_size >>> (
                  d_x_spin_buffer, d_y_spin_buffer, d_z_spin_buffer,
                  d_materials, d_heun_params,
                  d_x_spin_field, d_y_spin_field, d_z_spin_field,
                  d_x_external_field, d_y_external_field, d_z_external_field,
                  d_x_spin, d_y_spin, d_z_spin,
                  ::mp::dt, ::atoms::num_atoms
                  );

#ifdef CUDA_DEBUG
            std::cout << cu::atoms::x_spin_array[0] << " "
                      << cu::atoms::y_spin_array[0] << " "
                      << cu::atoms::z_spin_array[0] << std::endl;
#endif

            check_cuda_errors (__FILE__, __LINE__);

            cu::update_spin_fields ();

            check_cuda_errors (__FILE__, __LINE__);

            /*
             * TODO: Store delta s because of the renormalization
             */
            cu::llg::llg_heun_scheme <<< cu::grid_size, cu::block_size >>> (
                  d_x_spin, d_y_spin, d_z_spin,
                  d_materials, d_heun_params,
                  d_x_spin_field, d_y_spin_field, d_z_spin_field,
                  d_x_external_field, d_y_external_field, d_z_external_field,
                  d_x_spin_buffer, d_y_spin_buffer, d_z_spin_buffer,
                  ::mp::dt, ::atoms::num_atoms
                  );

            check_cuda_errors (__FILE__, __LINE__);

            /*
             * TODO: This copy can go away
             *       The buffer contains the old spin and spin contains
             *       the updated version.
             */
            thrust::copy (
                  cu::llg::x_spin_buffer_array.begin(),
                  cu::llg::x_spin_buffer_array.end(),
                  cu::atoms::x_spin_array.begin()
                  );
            thrust::copy (
                  cu::llg::y_spin_buffer_array.begin(),
                  cu::llg::y_spin_buffer_array.end(),
                  cu::atoms::y_spin_array.begin()
                  );
            thrust::copy (
                  cu::llg::z_spin_buffer_array.begin(),
                  cu::llg::z_spin_buffer_array.end(),
                  cu::atoms::z_spin_array.begin()
                  );
         }

         __global__ void llg_heun_step (
               double * x_spin, double * y_spin, double * z_spin,
               int * material_id,
               cu::heun_parameters_t * heun_parameters,
               double * x_sp_field, double * y_sp_field, double * z_sp_field,
               double * x_ext_field, double * y_ext_field, double * z_ext_field,
               double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
               double dt, size_t num_atoms
               )
         {

            for ( size_t atom = blockIdx.x * blockDim.x + threadIdx.x;
                  atom < num_atoms;
                  atom += blockDim.x * gridDim.x)
            {

               size_t mid = material_id[atom];

               double prefactor = heun_parameters[mid].prefactor;
               double lambdatpr = heun_parameters[mid].lambda_times_prefactor;

               double sx = x_spin[atom];
               double sy = y_spin[atom];
               double sz = z_spin[atom];

               double H_x = x_sp_field[atom] + x_ext_field[atom];
               double H_y = y_sp_field[atom] + y_ext_field[atom];
               double H_z = z_sp_field[atom] + z_ext_field[atom];

               //defining the s x h array and s x s x h :warray
               double sxh_x = sy * H_z - sz * H_y;
               double sxh_y = sz * H_x - sx * H_z;
               double sxh_z = sx * H_y - sy * H_x;

               //defining the sxsxh
               double sxsxh_x = sy * sxh_z - sz * sxh_y;
               double sxsxh_y = sz * sxh_x - sx * sxh_z;
               double sxsxh_z = sx * sxh_y - sy * sxh_x;

               //the saturation
               double mod_s = 0.0;

               //defining the Delta
               double Ds_x = - prefactor * sxh_x + lambdatpr * sxsxh_x;
               double Ds_y = - prefactor * sxh_y + lambdatpr * sxsxh_y;
               double Ds_z = - prefactor * sxh_z + lambdatpr * sxsxh_z;

               double new_spin_x = sx + Ds_x * dt;
               double new_spin_y = sy + Ds_y * dt;
               double new_spin_z = sz + Ds_z * dt;

               /*
                * TODO: Adjust delta s for the non norm conserving stuff
                */
               mod_s = 1.0f / sqrtf(
                     new_spin_x * new_spin_x +
                     new_spin_y * new_spin_y +
                     new_spin_z * new_spin_z);

               //normalized spins
               x_spin_prim[atom] = new_spin_x * mod_s;
               y_spin_prim[atom] = new_spin_y * mod_s;
               z_spin_prim[atom] = new_spin_z * mod_s;
            }

         }

         __global__ void llg_heun_scheme (
               double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
               int * material_id,
               cu::heun_parameters_t * heun_parameters,
               /*
                * receive spin init
                */
               double * x_sp_field, double * y_sp_field, double * z_sp_field,
               double * x_ext_field, double * y_ext_field, double * z_ext_field,
               /*
                * receive spin prima, read a write the final result to it
                */
               double * x_spin, double * y_spin, double * z_spin,
               double dt, size_t num_atoms
               )
         {

            for ( size_t atom = blockIdx.x * blockDim.x + threadIdx.x;
                  atom < num_atoms;
                  atom += blockDim.x * gridDim.x)
            {

               size_t mid = material_id[atom];

               double prefactor = heun_parameters[mid].prefactor;
               double lambdatpr = heun_parameters[mid].lambda_times_prefactor;

               //heun step array
               /*
                * TODO: Update this, read in the delta s
                */
               double Ds_x = (x_spin_prim[atom] - x_spin[atom]) / dt;
               double Ds_y = (y_spin_prim[atom] - y_spin[atom]) / dt;
               double Ds_z = (z_spin_prim[atom] - z_spin[atom]) / dt;

               //initial spins
               double spin_init_x = x_spin[atom];
               double spin_init_y = y_spin[atom];
               double spin_init_z = z_spin[atom];

               //update the spins
               double spin_x = x_spin_prim[atom];
               double spin_y = y_spin_prim[atom];
               double spin_z = z_spin_prim[atom];

               //the field
               double H_x = x_sp_field[atom] + x_ext_field[atom];
               double H_y = y_sp_field[atom] + y_ext_field[atom];
               double H_z = z_sp_field[atom] + z_ext_field[atom];

               //implementing the Heun's Scheme
               //cross product
               double SxH_x = spin_x * H_z - spin_z * H_y;
               double SxH_y = spin_z * H_x - spin_x * H_z;
               double SxH_z = spin_x * H_y - spin_y * H_x;

               double SxSxH_x = spin_y * SxH_z - spin_z * SxH_y;
               double SxSxH_y = spin_z * SxH_x - spin_x * SxH_z;
               double SxSxH_z = spin_x * SxH_y - spin_y * SxH_x;

               double DS_prime_x = - prefactor * SxH_x + lambdatpr * SxSxH_x;
               double DS_prime_y = - prefactor * SxH_y + lambdatpr * SxSxH_y;
               double DS_prime_z = - prefactor * SxH_z + lambdatpr * SxSxH_z;

               double S_x = spin_init_x + 0.5f * (Ds_x + DS_prime_x) * dt;
               double S_y = spin_init_y + 0.5f * (Ds_y + DS_prime_y) * dt;
               double S_z = spin_init_z + 0.5f * (Ds_z + DS_prime_z) * dt;

               float mods = 1.0f / sqrtf(
                     S_x*S_x + S_y*S_y + S_z*S_z);

               x_spin[atom] = mods * S_x;
               y_spin[atom] = mods * S_y;
               z_spin[atom] = mods * S_z;
            }
         }
      } /* llg */
   }

#endif

} // end of namespace cuda
