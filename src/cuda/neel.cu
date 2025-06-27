//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
// Reviewed: Andrea Meo 2022
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

//------------------------------------------------------------------------------
// Host function to calculate neel fields using gpu kernel
//------------------------------------------------------------------------------
void update_neel_fields()
{
   //Call kernel to calculate neel fields
   cu::update_neel_fields_kernel <<< cu::grid_size, cu::block_size >>> (
         cu::d_neel_tensor, cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin, cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field, ::atoms::num_atoms
      );

   check_cuda_errors(__FILE__, __LINE__);
}

__global__ void update_neel_fields_kernel (
   cu_real_t * d_neel_tensor,
   cu_real_t * x_spin,
   cu_real_t * y_spin,
   cu_real_t * z_spin,
   cu_real_t * x_sp_field,
   cu_real_t * y_sp_field,
   cu_real_t * z_sp_field,
   int n_atoms
   )
{
   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x) {

      const cu_real_t sx = x_spin[i];
      const cu_real_t sy = y_spin[i];
      const cu_real_t sz = z_spin[i];

      //Get index for atom inside tensor
      const int index = 9*i;

      // Second order neel
      cu_real_t hx = 2.0 * ( d_neel_tensor[index + 0] * sx +
                             d_neel_tensor[index + 1] * sy +
                             d_neel_tensor[index + 2] * sz );

      cu_real_t hy = 2.0 * ( d_neel_tensor[index + 3] * sx +
                             d_neel_tensor[index + 4] * sy +
                             d_neel_tensor[index + 5] * sz );

      cu_real_t hz = 2.0 * ( d_neel_tensor[index + 6] * sx +
                             d_neel_tensor[index + 7] * sy +
                             d_neel_tensor[index + 8] * sz );

      //Store net field
      x_sp_field[i] += hx;
      y_sp_field[i] += hy;
      z_sp_field[i] += hz;
   }

   return;
}

} // end of internal namespace

} // end of vcuda namespace

#endif
