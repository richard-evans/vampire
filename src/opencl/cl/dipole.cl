#include "cl_defs.h"
#include "material_type.h"

__kernel
void update_dipole_fields(const __global real_t *const restrict x_mag,
                          const __global real_t *const restrict y_mag,
                          const __global real_t *const restrict z_mag,
                          const __global real_t *const restrict x_coord,
                          const __global real_t *const restrict y_coord,
                          const __global real_t *const restrict z_coord,
                          const __global real_t *const restrict volume,
                          __global real_t *const restrict x_dip_field,
                          __global real_t *const restrict y_dip_field,
                          __global real_t *const restrict z_dip_field)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_CELLS; i+=gsz)
   {
      real_t mx = x_mag[i];
      real_t my = y_mag[i];
      real_t mz = z_mag[i];
      real_t cx = x_coord[i];
      real_t cy = y_coord[i];
      real_t cz = z_coord[i];

      real_t vol_prefac = - 4 * PI / (3 * volume[i]);
      real_t prefactor  = 1e23;

      real_t field_x = vol_prefac * mx;
      real_t field_y = vol_prefac * my;
      real_t field_z = vol_prefac * mz;

      for (size_t j=0; j<N_CELLS; ++j)
      {
         if (i==j) continue;

         real_t omx = x_mag[i];
         real_t omy = y_mag[i];
         real_t omz = z_mag[i];

         real_t dx = x_coord[j] - cx;
         real_t dy = y_coord[j] - cy;
         real_t dz = z_coord[j] - cz;

         real_t drij  = RSQRT(dx*dx + dy*dy + dz*dz);
         real_t drij3 = drij * drij * drij;

         real_t sdote = (omx * dx * drij +
                         omy * dy * drij +
                         omz * dz * drij);

         field_x += (3 * sdote * dx * drij - omx) * drij3;
         field_y += (3 * sdote * dy * drij - omy) * drij3;
         field_z += (3 * sdote * dz * drij - omz) * drij3;
      }

      x_dip_field[i] = prefactor * field_x;
      y_dip_field[i] = prefactor * field_y;
      z_dip_field[i] = prefactor * field_z;
   }
}

__kernel
void update_atm_dipole_fields(const __global real_t *const restrict x_cell_field,
                              const __global real_t *const restrict y_cell_field,
                              const __global real_t *const restrict z_cell_field,
                              __global real_t *const restrict x_dip_field,
                              __global real_t *const restrict y_dip_field,
                              __global real_t *const restrict z_dip_field,
                              const __global int *const restrict cell)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      int cid = cell[i];
      x_dip_field[i] = x_cell_field[cid];
      y_dip_field[i] = y_cell_field[cid];
      z_dip_field[i] = z_cell_field[cid];
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
void update_cell_magnetization(const __global real_t *const restrict x_spin,
                               const __global real_t *const restrict y_spin,
                               const __global real_t *const restrict z_spin,
                               const __global int *const restrict material,
                               const __global int *const restrict cell,
                               const __global material_parameters_t *const restrict material_params,
                               __global real_t *const restrict x_mag,
                               __global real_t *const restrict y_mag,
                               __global real_t *const restrict z_mag)
{
   size_t gsz = get_global_size(0);

   for (size_t i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      int mid = material[i];
      int cid = cell[i];
      real_t mu_s = material_params[mid].mu_s_si;

      atomic_add_global(&x_mag[cid], x_spin[i]*mu_s);
      atomic_add_global(&y_mag[cid], y_spin[i]*mu_s);
      atomic_add_global(&z_mag[cid], z_spin[i]*mu_s);
   }
}
