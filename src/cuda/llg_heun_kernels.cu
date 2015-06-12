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
	////////////////////////////////////////////
	// const gyro;
	// const alfa	
	//////////////////////////////////////////////
	int atom = blockIdx.x * blockDim.x + threadIdx.x;
	//the total field array
	float3 H	
	H_x[atom] = x_sp_field + x_ext_field;
	H_y[atom] = y_sp_field + y_ext_field;
	H_z[atom] = z_sp_field + z_ext_field; 

         //defining the s x h array and s x s x h :warray
	float3 sxh
	sxh.x = y_spin * H_z[atom] - z_spin * H_y[atom];
	sxh.y = z_spin * H_x[atom] - x_spin * H_z[atom];
	sxh.z = x_spin * H_y[atom] - y_spin * H_x[atom];
	
	//defining the sxsxh
	float3 sxsxh;	 
	sxsxh.x = y_spin * sxh.z - z_spin * sxh.y;
        sxsxh.y = z_spin * sxh.x - x_spin * sxh.z;
        sxsxh.z = x_spin * sxh.y - y_spin * sxh.x;
	
	//the saturation
	float mod_s = 0.0

	//defining the Delta
	float3 Ds
	Ds.x = gyro/(1 + alfa*alfa ) * (sxh.x + alfa*sxsxh.x);
        Ds.y = gyro/(1 + alfa*alfa ) * (sxh.y + alfa*sxsxh.y);
	Ds.z = gyro/(1 + alfa*alfa ) * (sxh.z + alfa*sxsxh.z);

	float3 new_spin
	new_spin.x = x_spin + Ds.x*dt; 
	new_spin.y = y_spin + Ds.y*dt;
	new_spin.z = z_spin + Ds.z*dt;

        mods = sqrtf(new_spin.x*new_spin.x + new_spin.y*new_spin.y + new_spin.z*new_spin.z);

	//normalized spins
	float3 m;
	m.x = new_spin.x/mods;
	m.y = new_spin.y/mods;
	m.z = new_spin.z/mods;
      }

   }
}
