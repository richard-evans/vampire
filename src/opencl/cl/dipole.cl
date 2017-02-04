#include "cl_defs.h"
#include "material_type.h"

#ifdef USE_VECTOR_TYPE
typedef real_t3 T;
#else
typedef real_t T;
#endif

__kernel
void update_dipole_fields(const __global T *const restrict mag,
                          const __global T *const restrict coord,
                          const __global real_t *const restrict volume,
                          __global T *const restrict dip_field)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_CELLS; i+=gsz)
   {
#ifdef USE_VECTOR_TYPE
      const real_t3 mi = mag[i];
      const real_t3 ci = coord[i];
#else
      const size_t xi = 3*i+0;
      const size_t yi = 3*i+1;
      const size_t zi = 3*i+2;

      const real_t3 mi = (real_t3)(mag[xi], mag[yi], mag[zi]);
      const real_t3 ci = (real_t3)(coord[xi], coord[yi], coord[zi]);
#endif

      const real_t vol_prefac = - 4 * PI / (3 * volume[i]);
      const real_t prefactor  = 1e23;

      real_t3 field = vol_prefac * mi;

      for (size_t j=0; j<N_CELLS; ++j)
      {
         if (i==j) continue;

#ifdef USE_VECTOR_TYPE
         const real_t3 mj = mag[j];
         const real_t3 cj = coord[j];
#else
         const size_t xj = 3*j+0;
         const size_t yj = 3*j+1;
         const size_t zj = 3*j+2;

         const real_t3 mj = (real_t3)(mag[xj], mag[yj], mag[zj]);
         const real_t3 cj = (real_t3)(coord[xj], coord[yj], coord[zj]);
#endif
         const real_t3 dX = cj - ci;

         const real_t drij  = RSQRT(dX.x*dX.x + dX.y*dX.y + dX.z*dX.z);
         const real_t drij3 = drij * drij * drij;

         const real_t3 sdote_vec = mj * dX * drij;
         const real_t sdote = sdote_vec.x + sdote_vec.y + sdote_vec.z;

         field += (3 * sdote * dX * drij - mj) * drij3;
      }

      field *= prefactor;

#ifdef USE_VECTOR_TYPE
      dip_field[i] = field;
#else
      dip_field[xi] = field.x;
      dip_field[yi] = field.y;
      dip_field[zi] = field.z;
#endif
   }
}

__kernel
void update_atm_dipole_fields(const __global T *const restrict cell_field,
                              __global T *const restrict dip_field,
                              const __global int *const restrict cell)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      const int cid = cell[i];

#ifdef USE_VECTOR_TYPE
      dip_field[i] = cell_field[i];
#else
      dip_field[3*i+0] = cell_field[3*cid+0];
      dip_field[3*i+1] = cell_field[3*cid+1];
      dip_field[3*i+2] = cell_field[3*cid+2];
#endif
   }
}

void atomic_add_global(volatile __global real_t *const source,
                       const real_t operand)
{
   union
   {
       uint_t i;
       real_t f;
   } newVal, prevVal;

   do
   {
      prevVal.f = *source;
      newVal.f = prevVal.f + operand;
   }
   while (ATOMIC_CMPXCHG((volatile __global uint *const)source, prevVal.i, newVal.i) != prevVal.i);
}

__kernel
void update_cell_magnetization(const __global T *const restrict spin,
                               const __global int *const restrict material,
                               const __global int *const restrict cell,
                               const __global material_parameters_t *const restrict material_params,
                               __global T *const restrict mag)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      const int mid = material[i];
      const int cid = cell[i];
      const real_t mu_s = material_params[mid].mu_s_si;

#ifdef USE_VECTOR_TYPE
      atomic_add_global(&mag[cid].x, spin[i].x*mu_s);
      atomic_add_global(&mag[cid].y, spin[i].y*mu_s);
      atomic_add_global(&mag[cid].z, spin[i].z*mu_s);
#else
      const size_t x = 3*i+0;
      const size_t y = 3*i+1;
      const size_t z = 3*i+2;

      atomic_add_global(&mag[3*cid+0], spin[x]*mu_s);
      atomic_add_global(&mag[3*cid+1], spin[y]*mu_s);
      atomic_add_global(&mag[3*cid+2], spin[z]*mu_s);
#endif
   }
}
