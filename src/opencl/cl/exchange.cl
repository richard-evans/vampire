#include "cl_defs.h"

#ifdef USE_VECTOR_TYPE
typedef real_t3 T;
#else
typedef real_t T;
#endif

__kernel
void calculate_exchange(const int exchange_type,
                        const __global real_t *const restrict Jxx,
                        const __global real_t *const restrict Jyy,
                        const __global real_t *const restrict Jzz,
                        const __global uint   *const restrict limits,
                        const __global uint   *const restrict neighbours,
                        const __global T *const restrict spin,
                        __global T *const restrict total_spin_field)
{
   size_t gsz = get_global_size(0);

   switch (exchange_type)
   {
   case 0:
      // Isotropic
      for (size_t i=get_global_id(0); i<N; i+=gsz)
      {
         real_t3 sum = (real_t3)(0.0, 0.0, 0.0);

#ifdef USE_VECTOR_TYPE
         for (uint j=limits[i]; j<limits[i+1]; ++j)
         {
            sum += Jxx[j] * spin[neighbours[j]];
         }
         total_spin_field[i] = sum;
#else
         const size_t xi = 3*i+0;
         const size_t yi = 3*i+1;
         const size_t zi = 3*i+2;

         for (uint j=limits[i]; j<limits[i+1]; ++j)
         {
            size_t nx = 3*neighbours[j]+0;
            size_t ny = 3*neighbours[j]+1;
            size_t nz = 3*neighbours[j]+2;
            real_t3 snj = (real_t3)(spin[nx], spin[ny], spin[nz]);
            sum += Jxx[j] * snj;
         }

         total_spin_field[xi] = sum.x;
         total_spin_field[yi] = sum.y;
         total_spin_field[zi] = sum.z;
#endif // USE_VECTOR_TYPE
      }
      break;

   case 1:
      // Vector
      for (unsigned i=0; i<N; i+=gsz)
      {
         real_t3 sum = (real_t3)(0.0, 0.0, 0.0);

#ifdef USE_VECTOR_TYPE
         for (uint j=limits[i]; j<limits[i+1]; ++j)
         {
            real_t3 J = (real_t3)(Jxx[j], Jyy[j], Jzz[j]);
            sum += J * spin[neighbours[j]];
         }
         total_spin_field[i] = sum;
#else
         const size_t xi = 3*i+0;
         const size_t yi = 3*i+1;
         const size_t zi = 3*i+2;

         for (uint j=limits[i]; j<limits[i+1]; ++j)
         {
            size_t nx = 3*neighbours[j]+0;
            size_t ny = 3*neighbours[j]+1;
            size_t nz = 3*neighbours[j]+2;
            real_t3 snj = (real_t3)(spin[nx], spin[ny], spin[nz]);
            real_t3 J = (real_t3)(Jxx[j], Jyy[j], Jzz[j]);
            sum += J * snj;
         }
         total_spin_field[x] = sum.x;
         total_spin_field[y] = sum.y;
         total_spin_field[z] = sum.z;
#endif // USE_VECTOR_TYPE
      }
   }
}
