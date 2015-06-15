#include "internal.hpp"
#include <cuda.h>
namespace cuda {
   namespace internal {

      __global__ void llg_heun_first_kernel (
            double * x_spin, double * y_spin, double * z_spin,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            double dt
            )
      {
         int atom = blockIdx.x * blockDim.x + threadIdx.x;
         //the total field array

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
         Ds.x = gyro/(1 + alfa*alfa ) * (sxh.x + alfa*sxsxh.x);
         Ds.y = gyro/(1 + alfa*alfa ) * (sxh.y + alfa*sxsxh.y);
         Ds.z = gyro/(1 + alfa*alfa ) * (sxh.z + alfa*sxsxh.z);

         float3 new_spin;
         new_spin.x = spin.x + Ds.x*dt;
         new_spin.y = spin.y + Ds.y*dt;
         new_spin.z = spin.z + Ds.z*dt;

         mods = sqrtf(new_spin.x*new_spin.x + new_spin.y*new_spin.y + new_spin.z*new_spin.z);

         //normalized spins
         float3 m;
         m.x = new_spin.x/mods;
         m.y = new_spin.y/mods;
         m.z = new_spin.z/mods;
      }

   }
}
