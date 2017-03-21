//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "cl_defs.h"

// performs matrix multiplication using a CSR format matrix for the exchange calculation
// if USE_VECTOR_TYPE is defined, A should be a buffer of 3 component vectors, for Jxx,
// Jyy and Jzz. Then R.x = matmul(A.x, V.x) and similarly for y and z.

__kernel
void matmul(const __global real_t *const restrict A,  /* matrix non-zero values  */
            const __global uint   *const restrict IA, /* matrix row ptr */
            const __global uint   *const restrict JA, /* col idx */
            const __global real_t *const restrict V,  /* vector to be multiplied */
                  __global real_t *const restrict R)  /* resultant vector */
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<NUM_ATOMS; i+=gsz)
   {
      real_t3 sum = (real_t3)(0.0, 0.0, 0.0);

      for (uint j=IA[i]; j<IA[i+1]; ++j)
      {
#if EXCH_TYPE == 0
         sum += A[j] * vload3(JA[j], V);

#elif EXCH_TYPE == 1
         sum += vload3(j, A) * vload3(JA[j], V);

#else
#error EXCH_TYPE not supported
#endif
      }

      vstore3(sum, i, R);
   }
}
