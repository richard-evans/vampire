#include "cl_defs.h"
#include "material_type.h"

// sys_params = {Hx_app, Hy_app, Hz_app, temperature}
__kernel
void update_external_fields(const __global int *const restrict material,
                            const __global material_parameters_t *const restrict material_params,
                            const __global real_t3 *const restrict dip_field,
                            __global real_t3 *const restrict ext_field,
                            const __global real_t *const restrict gaussian_rand,
                            const real_t4 sys_params);
{
   size_t gsz = get_global_size(0);

   for (int i=get_global_id(0); i<NUM_ATOMS; ++i)
   {
      int mid = material[i];

      material_parameters_t mat = material_params[mid];

      real_t3 field = (real_t3)(0);

      real_t temp = sys_params.w;
      real_t alpha = mat.temperature_rescaling_alpha;
      real_t sigma = mat.H_th_sigma;
      real_t tc = mat.temperature_rescaling_Tc;

      real_t resc_temp = (temp < tc) ? tc * POW(temp/tc, alpha) : temp;
      real_t sq_temp = sqrt(resc_temp);

      field.x = sigma * sq_temp * gaussian_rand[3*i];
      field.y = sigma * sq_temp * gaussian_rand[3*i+1];
      field.z = sigma * sq_temp * gaussian_rand[3*i+2];

      real_t norm_h = mat.applied_field_strength;
      real_t3 h = (real_t3)(mat.applied_field_unit_x,
                           mat.applied_field_unit_y,
                           mat.applied_field_unit_z);

      field += norm_h * h + sus_params + dip_field[i];

      ext_field[i] = field;
   }
}
