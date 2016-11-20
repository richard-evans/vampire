#include "material_type.h"

#ifdef OPENCL_DP
typedef double real_t
#else
typedef float  real_t
#endif

#ifdef OPENCL_USE_NATIVE_FUNCTIONS
#define POW(x, y) native_powr(x, y)
#else
#define POW(x, y) pow(x, y)
#endif

__kernel
void update_external_fields(const __global int *material,
                            const __global material_params_t material_params,
                            const __global real_t *x_dip_field,
                            const __global real_t *y_dip_field,
                            const __global real_t *z_dip_field,
                            __global real_t *x_ext_field,
                            __global real_t *y_ext_field,
                            __global real_t *z_ext_field,
                            const __global real_t *gaussian_rand,
                            __constant real_t global_temperature,
                            __constant real_t Hx_app,
                            __constant real_t Hy_app,
                            __constant real_t Hz_app)
{
   size_t gsz = get_global_size(0);

   for (int i=get_global_id(0); i<N_ATOMS; ++i)
   {
      int mid = material[i];

      material_params_t mat material_params[mid];

      real_t field_x = 0;
      real_t field_y = 0;
      real_t field_z = 0;

      real_t temp = global_termperature;
      real_t alpha = mat.temperature_rescaling_alpha;
      real_t sigma = mat.H_th_sigma;
      real_t tc = mat.temperature_rescaling_Tc;

      real_t resc_temp = (temp < tc) ? tc * POW(temp/tc, alpha) : temp;
      real_t sq_temp = sqrt(resc_temp);

      field_x = sigma * sq_temp * gaussian_rand[3*i];
      field_y = sigma * sq_temp * gaussian_rand[3*i+1];
      field_z = sigma * sq_temp * gaussian_rand[3*i+2];

      real_t norm_h = mat.applied_field_strength;
      real_t hx = mat.applied_field_unit_x;
      real_t hy = mat.applied_field_unit_y;
      real_t hz = mat.applied_field_unit_z;

      field_x += norm_h * hx + Hx_app + x_dip_field[i];
      field_y += norm_h * hy + Hy_app + y_dip_field[i];
      field_z += norm_h * hz + Hz_app + z_dip_field[i];

      x_ext_field[i] = field_x;
      y_ext_field[i] = field_y;
      z_ext_field[i] = field_z;
   }
}
