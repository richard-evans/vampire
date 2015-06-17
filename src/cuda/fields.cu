#include "internal.hpp"

namespace cuda
{
   namespace internal
   {

      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material, size_t * cell,
            double * x_sp_field, double * y_sp_field, double * z_sp_field
            )
      {
         size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      }

   } /* internal */
} /* cuda */

