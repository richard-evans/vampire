#include "cl_defs.h"

__kernel
void matmul(const __global real_t *Mvals,
            const __global uint   *Mrow_ptr,
            const __global uint   *Mcol_idx,
            const __global real_t *V,
            __global   real_t *Ret)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N; ++i)
   {
      real_t sum = 0;

      for (uint j=Mrow_ptr[i]; j<Mrow_ptr[i+1]; ++j) {
         sum += Mvals[j] * V[Mcol_idx[j]];
      }

      Ret[i] = sum;
   }
}
