#include <cuda.h>
#include "internal.hpp"
namespace cuda {
 	namespace internal {

	__global__ void llg_heun_scheme( double * x_spin, double * y_spin, double * z_spin,
				         double * x_sp_field, double * y_sp_field, double * z_sp_field,
					 double * x_ext_field, double * y_ext_field, double * z_ext_field,   
					 double * x_new_spin, double * y_new_spin, double * z_new_spin,
					 double * Ds.x, double * Ds.y, double * Ds.z, double dt){

	int atom = blockIdx.x * blockDim.x + threadIdx.x; 
	
	//heun step array
	float3 Ds;
	Ds.x = Ds.x[atom];
	Ds.y = Ds.y[atom];
	Ds.z = Ds.z[atom];

	//initial spins
	float3 spin_init;
	spin_init = x_spin[atom];
	spin_init = y_spin[atom];
	spin_init = z_spin[atom];

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
	SxH.y = spin.z * H.x - spin.x * H.z
	SxH.z = spin.x * H.y - spin.y * H.x;
	
	float3 SxSxH;
	SxSxH.x = spin.y * SxH.z - spin.z * SxH.y;
	SxSxH.y = spin.z * SxH.x - spin.x * SxH.z;
	SxSxH.z = spin.x * SxH.y - spin.y * SxH.x;
	
	float3 DS_prime;
	DS_prime.x = -gyro/(1 + alfa*alfa ) * (SxH.x + alfa*SxSxH.x);
	DS_prime.y = -gyro/(1 + alfa*alfa ) * (SxH.y + alfa*SxSxH.y);
	DS_prime.z = -gyro/(1 + alfa*alfa ) * (SxH.z + alfa*SxSxH.z);
	
	float3 S;
	S.x = spin.x + 0.5 * (Ds.x + DS_prime.x) * dt;
	S.y = spin.y + 0.5 * (Ds.y + DS_prime.y) * dt;
	S.z = spin.z + 0.5 * (Ds.z + DS_prime.z) * dt;
	
	float mods =0.0;
	mods = 1/sqrtf(S.x*S.x + S.y*S.y + S.z*S.z);
	
	float3 Spin;
	Spin.x = mods * S.x;
	Spin.y = mods * S.y;
	Spin.z = mods * S.z;
	}    
}   
  }
