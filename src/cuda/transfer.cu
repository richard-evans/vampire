//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//------------------------------------------------------------------------------
//
// C++ standard library headers

// Vampire headers
#include "cuda.hpp"


// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "dipole.hpp"
#include "internal.hpp"
#include "typedefs.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

//------------------------------------------------------------------------------
// Wrapper function to transfer spin data from GPU to CPU
//------------------------------------------------------------------------------
void transfer_spin_positions_from_gpu_to_cpu(){

   cudaMemcpy(::atoms::x_spin_array.data(), cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);
   cudaMemcpy(::atoms::y_spin_array.data(), cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);
   cudaMemcpy(::atoms::z_spin_array.data(), cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyDeviceToHost);
   
   /*
   thrust::copy(internal::atoms::x_spin_array.begin(),internal::atoms::x_spin_array.end(),::atoms::x_spin_array.begin());
   thrust::copy(internal::atoms::y_spin_array.begin(),internal::atoms::y_spin_array.end(),::atoms::y_spin_array.begin());
   thrust::copy(internal::atoms::z_spin_array.begin(),internal::atoms::z_spin_array.end(),::atoms::z_spin_array.begin());
   */
   return;

}

//------------------------------------------------------------------------------
// Wrapper function to transfer dipole field data from CPU to GPU
//------------------------------------------------------------------------------
void transfer_dipole_fields_from_cpu_to_gpu(){

   cudaMemcpy(::atoms::x_spin_array.data(), cu::atoms::d_x_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyHostToDevice);
   cudaMemcpy(::atoms::y_spin_array.data(), cu::atoms::d_y_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyHostToDevice);
   cudaMemcpy(::atoms::z_spin_array.data(), cu::atoms::d_z_spin, ::atoms::num_atoms * sizeof(cu::cu_real_t), cudaMemcpyHostToDevice);

   /*
   thrust::copy(::dipole::atom_dipolar_field_array_x.begin(),::dipole::atom_dipolar_field_array_x.end(), cu::x_dipolar_field_array.begin());
   thrust::copy(::dipole::atom_dipolar_field_array_y.begin(),::dipole::atom_dipolar_field_array_y.end(), cu::y_dipolar_field_array.begin());
   thrust::copy(::dipole::atom_dipolar_field_array_z.begin(),::dipole::atom_dipolar_field_array_z.end(), cu::z_dipolar_field_array.begin());
   */
   return;

}

} // end of vcuda namespace

#endif
