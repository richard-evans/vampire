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
#include "internal.hpp"
#include "exchange_fields.hpp"

#ifdef CUDA
namespace cu = ::vcuda::internal;
#endif

namespace vcuda{

   //--------------------------------------------------------------------------
   // Function to perform a single heun step
   //--------------------------------------------------------------------------
   void llg_heun(){

      cu::exchange::calculate_exchange_fields();


#ifdef CUDA
      /* set up and call the kernels */
      /* assume that you have the data already
       * in the device */

#endif

      return;
   }

#ifdef CUDA

   namespace internal {

      __global__ void llg_heun_first_kernel (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material_id,
            cu::heun_parameters_t * heun_parameters,
            double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
            double * x_delta_spin, double * y_delta_spin, double * z_delta_spin,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
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

            float3 spin;

            spin.x = x_spin[atom];
            spin.y = y_spin[atom];
            spin.z = z_spin[atom];

            float3 H;
            H.x = x_sp_field[atom] + x_ext_field[atom];
            H.y = y_sp_field[atom] + y_ext_field[atom];
            H.z = z_sp_field[atom] + z_ext_field[atom];

            //defining the s x h array and s x s x h :warray
            float3 sxh;
            sxh.x = spin.y * H.z - spin.z * H.y;
            sxh.y = spin.z * H.x - spin.x * H.z;
            sxh.z = spin.x * H.y - spin.y * H.x;

            //defining the sxsxh
            float3 sxsxh;
            sxsxh.x = spin.y * sxh.z - spin.z * sxh.y;
            sxsxh.y = spin.z * sxh.x - spin.x * sxh.z;
            sxsxh.z = spin.x * sxh.y - spin.y * sxh.x;

            //the saturation
            float mod_s = 0.0;

            //defining the Delta
            float3 Ds;
            Ds.x = -prefactor * sxh.x + lambdatpr * sxsxh.x;
            Ds.y = -prefactor * sxh.y + lambdatpr * sxsxh.y;
            Ds.z = -prefactor * sxh.z + lambdatpr * sxsxh.z;

            float3 new_spin;
            new_spin.x = spin.x + Ds.x*dt;
            new_spin.y = spin.y + Ds.y*dt;
            new_spin.z = spin.z + Ds.z*dt;

            mod_s = 1.0f/sqrtf(new_spin.x*new_spin.x + new_spin.y*new_spin.y + new_spin.z*new_spin.z);

            //normalized spins
            float3 m;
            m.x = new_spin.x*mod_s;
            m.y = new_spin.y*mod_s;
            m.z = new_spin.z*mod_s;

            // Update thetemporary buffers
         }

      }

      __global__ void llg_heun_scheme(
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material_id,
            cu::heun_parameters_t * heun_parameters,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            double * x_spin_prim, double * y_spin_prim, double * z_spin_prim,
            double * x_delta_spin, double * y_delta_spin, double * z_delta_spin,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            double * x_new_spin, double * y_new_spin, double * z_new_spin,
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
            float3 Ds;
            Ds.x = x_delta_spin[atom];
            Ds.y = y_delta_spin[atom];
            Ds.z = z_delta_spin[atom];

            //initial spins
            float3 spin_init;
            spin_init.x = x_spin[atom];
            spin_init.y = y_spin[atom];
            spin_init.z = z_spin[atom];

            //update the spins
            float3 spin;
            spin.x = x_new_spin[atom];
            spin.y = y_new_spin[atom];
            spin.z = z_new_spin[atom];

            //the field
            float3 H;
            H.x = x_sp_field[atom] + x_ext_field[atom];
            H.y = y_sp_field[atom] + y_ext_field[atom];
            H.z = z_sp_field[atom] + z_ext_field[atom];

            //implementing the Heun's Scheme
            //cross product
            float3 SxH;
            SxH.x = spin.x * H.z - spin.z * H.y;
            SxH.y = spin.z * H.x - spin.x * H.z;
            SxH.z = spin.x * H.y - spin.y * H.x;

            float3 SxSxH;
            SxSxH.x = spin.y * SxH.z - spin.z * SxH.y;
            SxSxH.y = spin.z * SxH.x - spin.x * SxH.z;
            SxSxH.z = spin.x * SxH.y - spin.y * SxH.x;

            float3 DS_prime;
            DS_prime.x = -prefactor * SxH.x + lambdatpr * SxSxH.x;
            DS_prime.y = -prefactor * SxH.y + lambdatpr * SxSxH.y;
            DS_prime.z = -prefactor * SxH.z + lambdatpr * SxSxH.z;

            float3 S;
            S.x = spin_init.x + 0.5f * (Ds.x + DS_prime.x) * dt;
            S.y = spin_init.y + 0.5f * (Ds.y + DS_prime.y) * dt;
            S.z = spin_init.z + 0.5f * (Ds.z + DS_prime.z) * dt;

            float mods = 1.0f/sqrtf(S.x*S.x + S.y*S.y + S.z*S.z);

            float3 Spin;
            Spin.x = mods * S.x;
            Spin.y = mods * S.y;
            Spin.z = mods * S.z;
         }
      }
   }

#endif

} // end of namespace cuda

