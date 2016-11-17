// Performs DIA format sparse matrix dense vector multiplication
// Mv = u -> matmul(u, M, v)

#ifdef OPENCL_DP
typedef double real_t;
#else
typedef float  real_t;
#endif

__kernel
void matmul(__global real_t *ret,
            const __global real_t *mvals,
            const __global real_t *vec)
{
   const size_t gsz = get_global_size(0);

   const uint offset = NDIAGS - (NDIAGS/2 + 1);

   for (uint i=get_global_id(0); i<N; i+=gsz)
   {
      real_t sum = 0;

      for (uint j=0; j<NDIAGS; ++j)
      {
         long vidx = i + j - offset;
         real_t velem = (0=<vidx && vidx<N) ? vec[vidx] : 0;
         sum += mvals[i*NDIAGS+j] * velem;
      }

      ret[i] = sum;
   }
}
