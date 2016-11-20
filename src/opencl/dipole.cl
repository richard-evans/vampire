#ifdef OPENCL_DP
typedef double real_t;
typedef ulong  uint_t;
#else
typedef float  real_t;
typedef uint   uint_t;
#endif

#ifdef OPENCL_USE_NATIVE_FUNCTIONS
#define RSQRT(x) native_rsqrt(x)
#else
#define RSQRT(c) rsqrt(x)
#endif

#define PI 3.14159265358979323846

__kernel
void update_dipole_fields(const __global real_t *x_mag,
                          const __global real_t *y_mag,
                          const __global real_t *z_mag,
                          const __global real_t *x_coord,
                          const __global real_t *y_coord,
                          const __global real_t *z_coord,
                          const __global real_t *volume,
                          __global real_t *x_dip_field,
                          __global real_t *y_dip_field,
                          __global real_t *z_dip_field)
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

__kernel update_atm_dipole_fields(const __global real_t *x_cell_field,
                                  const __global real_t *y_cell_field,
                                  const __global real_t *z_cell_field,
                                  __global real_t *x_dip_field,
                                  __global real_t *y_dip_field,
                                  __global real_t *z_dip_field,
                                  const __global int *cell)
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

void atomic_add_global(volatile __global real_t *source, const real_t operand)
{
   union
   {
      uint_t intVal;
      real_t floatVal;
   } newVal;
   union
   {
      uint_t intVal;
      real_t floatVal;
   } prevVal;

   do
   {
      prevVal.floatVal = *source;
      newVal.floatVal = prevVal.floatVal + operand;
   }
   while (atomic_cmpxchg((volatile __global uint *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

__kernel
void update_cell_magnetization(const __global real_t *x_spin,
                               const __global real_t *y_spin,
                               const __global real_t *z_spin,
                               const __global int *material,
                               const __global int *cell,
                               const __global material_params_t *material_params,
                               __global real_t *x_mag,
                               __global real_t *y_mag,
                               __global real_t *z_mag)
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
