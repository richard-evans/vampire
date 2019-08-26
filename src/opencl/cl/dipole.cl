//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "cl_defs.h"
#include "material_type.h"

__kernel
void update_dipole_fields(const __global vec_t  *const restrict mag,
                          const __global vec_t  *const restrict coord,
                          const __global real_t *const restrict volume,
                                __global vec_t  *const restrict dip_field)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_CELLS; i+=gsz)
   {
      const real_t3 mi = VEC_LOAD(mag, i);
      const real_t3 ci = VEC_LOAD(coord, i);

      const real_t vol_prefac = - 4 * PI / (3 * volume[i]);
      const real_t prefactor  = 1e23;

      real_t3 field = vol_prefac * mi;

      for (size_t j=0; j<N_CELLS; ++j)
      {
         if (i==j) continue;

         const real_t3 mj = VEC_LOAD(mag, j);
         const real_t3 cj = VEC_LOAD(coord, j);

         const real_t3 dX = cj - ci;

         const real_t drij  = RSQRT(dot(dX, dX));
         const real_t drij3 = drij * drij * drij;

         const real_t3 sdote_vec = mj * dX * drij;
         const real_t sdote = sdote_vec.x + sdote_vec.y + sdote_vec.z;

         field += (3 * sdote * dX * drij - mj) * drij3;
      }

      field *= prefactor;

      VEC_STORE(dip_field, i, field);
   }
}

__kernel
void update_atm_dipole_fields(const __global vec_t *const restrict cell_field,
                                    __global vec_t *const restrict dip_field,
                              const __global int   *const restrict cell)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      const int cid = cell[i];
      VEC_STORE(dip_field, i, VEC_LOAD(cell_field, cid));
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
void update_cell_magnetization(const __global vec_t *const restrict spin,
                               const __global int   *const restrict material,
                               const __global int   *const restrict cell,
                               const __global material_parameters_t *const restrict material_params,
                                     __global vec_t *const restrict mag)
{
   const size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      const int mid = material[i];
      const int cid = cell[i];
      const real_t mu_s = material_params[mid].mu_s_si;

#ifdef OPENCL_USE_VECTOR_TYPE
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
