#include "cl_defs.h"
#include "material_type.h"

__kernel
void update_nexch_spin_fields(const __global int *material,
                              const __global material_parameters_t *material_params,
                              const __global real_t *x_spin,
                              const __global real_t *y_spin,
                              const __global real_t *z_spin,
                              __global real_t *x_sp_field,
                              __global real_t *y_sp_field,
                              __global real_t *z_sp_field)
{
   size_t gsz = get_global_size(0);

   for (uint i=get_global_id(0); i<N_ATOMS; i+=gsz)
   {
      int mid = material[i];

      material_parameters_t mat = material_params[mid];

      real_t field_x = 0;
      real_t field_y = 0;
      real_t field_z = 0;

      real_t sx = x_spin[i];
      real_t sy = y_spin[i];
      real_t sz = z_spin[i];

      real_t ku = material.ku;
      field_z = - 2 * ku * sz;

      real_t ex = material.anisotropy_unit_x;
      real_t ey = material.anisotropy_unit_y;
      real_t ez = material.anisotropy_unit_z;

      real_t sdote  = sx*ex + sy*ey + sz*ez;
      real_t sdote3 = sdote * sdote * sdote;
      real_t sdote5 = sdote3 * sdote * sdote;

      real_t scale = real_t(2) / 3;

      real_t k2 = material.sh2;
      real_t k4 = material.sh4;
      real_t k6 = material.sh6;

      real_t ek2 = k2 * 3 * sdote;
      real_t ek4 = k4 * 0.125 * (140 * sdote3 - 60 * sdote);
      real_t ek6 = k6 * 0.0625 * (1386*sdote5 - 1260*sdote3 + 210*sdote);

      field_x += scale * ex * (ek2 + ek4 + ek6);
      field_y += scale * ey * (ek2 + ek4 + ek6);
      field_z += scale * ez * (ek2 + ek4 + ek6);

      x_sp_field[i] = field_x;
      y_sp_field[i] = field_y;
      z_sp_field[i] = field_z;
   }
}
