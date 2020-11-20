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
void update_nexch_spin_fields(const __global int *const restrict material,
                              const __global material_parameters_t *const restrict material_params,
                              const __global vec_t *const restrict spin,
                                    __global vec_t *const restrict sp_field)
{
   const size_t gsz = get_global_size(0);

   for (uint i=get_global_id(0); i<NUM_ATOMS; i+=gsz)
   {
      const real_t3 S = VEC_LOAD(spin, i);

      const int mid = material[i];

      const material_parameters_t mat = material_params[mid];

      real_t3 field = (real_t3)(0.0, 0.0, -2.0*mat.ku*S.z);

      const real_t3 e = (real_t3)(mat.anisotropy_unit_x,
                                  mat.anisotropy_unit_y,
                                  mat.anisotropy_unit_z);

      const real_t sdote = dot(S, e);
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

      VEC_INCR(sp_field, i, field);
   }
}
