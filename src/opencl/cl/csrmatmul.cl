//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "cl_defs.h"

#ifdef USE_VECTOR_TYPE
typedef real_t3 T;
#else
typedef real_t T;
#endif

// performs matrix multiplication using a CSR format matrix for the exchange calculation
// if USE_VECTOR_TYPE is defined, A should be a buffer of 3 component vectors, for Jxx,
// Jyy and Jzz. Then R.x = matmul(A.x, V.x) and similarly for y and z.

__kernel
void matmul(const __global T    *const restrict A,  /* matrix non-zero values  */
            const __global uint *const restrict IA, /* matrix row ptr */
            const __global uint *const restrict JA, /* col idx */
            const __global T    *const restrict V,  /* vector to be multiplied */
                  __global T    *const restrict R,  /* resultant vector */
            const uint N)                           /* length of vector */
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N; i+=gsz)
   {
      real_t3 sum = (real_t3)(0.0, 0.0, 0.0);

#ifdef USE_VECTOR_TYPE
      for (uint j=IA[i]; j<IA[i+1]; ++j)
      {
         sum += A[j] * V[JA[j]];
      }
      R[i] = sum;
#else
      const size_t xi = 3*i+0;
      const size_t yi = 3*i+1;
      const size_t zi = 3*i+2;

      for (uint j=IA[i]; j<IA[i+1]; ++j)
      {
         const size_t nx = 3*JA[j]+0;
         const size_t ny = 3*JA[j]+1;
         const size_t nz = 3*JA[j]+2;
         const real_t3 vnj = (real_t3)(V[nx], V[ny], V[nz]);
         sum += A[j] * vnj;
      }

      R[xi] = sum.x;
      R[yi] = sum.y;
      R[zi] = sum.z;
#endif // USE_VECTOR_TYPE
   }
}
