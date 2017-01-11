#include "cl_defs.h"
#include "material_type.h"

// sys_params = {Hx_app, Hy_app, Hz_app, temperature}
__kernel
void update_external_fields(const __global int *const restrict material,
                            const __global material_parameters_t *const restrict material_params,
                            const __global real_t *const restrict dip_field,
                            __global real_t *const restrict ext_field,
                            const __global real_t *const restrict gaussian_rand,
                            const real_t4 sys_params)
{
   size_t gsz = get_global_size(0);

   for (int i=get_global_id(0); i<NUM_ATOMS; ++i)
   {
      const size_t x = 3*i+0;
      const size_t y = 3*i+1;
      const size_t z = 3*i+2;

      const int mid = material[i];

      const material_parameters_t mat = material_params[mid];

      const real_t3 grands = (real_t3)(gaussian_rand[x], gaussian_rand[y], gaussian_rand[z]);

      const real_t temp = sys_params.w;
      const real_t alpha = mat.temperature_rescaling_alpha;
      const real_t sigma = mat.H_th_sigma;
      const real_t tc = mat.temperature_rescaling_Tc;

      const real_t resc_temp = (temp < tc) ? tc * POW(temp/tc, alpha) : temp;
      const real_t sq_temp = sqrt(resc_temp);

      real_t3 field = sigma * sq_temp * grands;

      const real_t norm_h = mat.applied_field_strength;
      const real_t3 h = (real_t3)(mat.applied_field_unit_x,
                                  mat.applied_field_unit_y,
                                  mat.applied_field_unit_z);

      field += norm_h * h + sys_params.xyz;

      ext_field[x] = field.x + dip_field[x];
      ext_field[y] = field.y + dip_field[y];
      ext_field[z] = field.z + dip_field[z];
   }
}
