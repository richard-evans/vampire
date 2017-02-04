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

#ifdef USE_VECTOR_TYPE
typedef real_t3 T;
#else
typedef real_t T;
#endif

__kernel
void update_nexch_spin_fields(const __global int *const restrict material,
                              const __global material_parameters_t *const restrict material_params,
                              const __global T *const restrict spin,
                              __global T *const restrict sp_field)
{
   const size_t gsz = get_global_size(0);

   for (uint i=get_global_id(0); i<NUM_ATOMS; i+=gsz)
   {
#ifdef USE_VECTOR_TYPE
      const real_t3 S = spin[i];
#else
      const size_t x = 3*i+0;
      const size_t y = 3*i+1;
      const size_t z = 3*i+2;

      const real_t3 S = (real_t3)(spin[x], spin[y], spin[z]);
#endif

      const int mid = material[i];

      const material_parameters_t mat = material_params[mid];

      real_t3 field = (real_t3)(0.0, 0.0, -2.0*mat.ku*S.z);

      const real_t3 e = (real_t3)(mat.anisotropy_unit_x,
                                  mat.anisotropy_unit_y,
                                  mat.anisotropy_unit_z);

      const real_t3 tmp = S * e;
      const real_t sdote  = tmp.x + tmp.y + tmp.z;
      const real_t sdote3 = sdote * sdote * sdote;
      const real_t sdote5 = sdote3 * sdote * sdote;

      const real_t scale = 2.0 / 3.0;

      const real_t k2 = mat.sh2;
      const real_t k4 = mat.sh4;
      const real_t k6 = mat.sh6;

      const real_t ek2 = k2 * 3 * sdote;
      const real_t ek4 = k4 * 0.125 * (140 * sdote3 - 60 * sdote);
      const real_t ek6 = k6 * 0.0625 * (1386*sdote5 - 1260*sdote3 + 210*sdote);

      const real_t ek_sum = ek2 + ek4 + ek6;
      field += scale * e * ek_sum;

#ifdef USE_VECTOR_TYPE
      sp_field[i] = field;
#else
      sp_field[x] = field.x;
      sp_field[y] = field.y;
      sp_field[z] = field.z;
#endif
   }
}
