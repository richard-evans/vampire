#include "cl_defs.h"
#include "material_type.h"

__kernel
void update_dipole_fields(const __global real_t3 *const restrict mag,
                          const __global real_t3 *const restrict coord,
                          const __global real_t *const restrict volume,
                          __global real_t3 *const restrict dip_field)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_CELLS; i+=gsz)
   {
      real_t3 m = mag[i];
      real_t3 c = coord[i];

      real_t vol_prefac = - 4 * PI / (3 * volume[i]);
      real_t prefactor  = 1e23;

      real_t3 field = vol_prefac * m;

      for (size_t j=0; j<N_CELLS; ++j)
      {
         if (i==j) continue;

         real_t3 om = mag[i];

         real_t3 dX = coord[j] - x;

         real_t drij  = RSQRT(dX.x*dX.x + dX.y*dX.y + dX.z*dX.z);
         real_t drij3 = drij * drij * drij;

         real_t3 sdote_vec = om * dx * drij;
         real_t sdote = sdote_vec.x + sdote_vec.y + sdote_vec.z;

         field += (3 * sdote * dX * drij - om) * drij3;
      }

      dip_field[i] = prefactor * field;
   }
}

__kernel
void update_atm_dipole_fields(const __global real_t3 *const restrict cell_field,
                              __global real_t3 *const restrict dip_field,
                              const __global int *const restrict cell)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      int cid = cell[i];
      dip_filed[i] = cell_field[cid];
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
void update_cell_magnetization(const __global real_t3 *const restrict spin,
                               const __global int *const restrict material,
                               const __global int *const restrict cell,
                               const __global material_parameters_t *const restrict material_params,
                               __global real_t3 *const restrict mag)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      int mid = material[i];
      int cid = cell[i];
      real_t mu_s = material_params[mid].mu_s_si;

      atomic_add_global(&mag[cid].x, spin[i].x*mu_s);
      atomic_add_global(&mag[cid].y, spin[i].y*mu_s);
      atomic_add_global(&mag[cid].z, spin[i].z*mu_s);
   }
}
