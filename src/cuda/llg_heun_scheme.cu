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
	
	
	}    
}   
  }
