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
void update_external_fields(const __global int *const restrict material,
                            const __constant material_parameters_t *const restrict material_params,
                            const __global T *const restrict dip_field,
                            __global T *const restrict ext_field,
                            const __global real_t *const restrict gaussian_rand,
                            const real_t Hx,
                            const real_t Hy,
                            const real_t Hz,
                            const real_t temp)
{
   size_t gsz = get_global_size(0);

   const real_t3 Happ = (real_t3)(Hx, Hy, Hz);

   for (int i=get_global_id(0); i<NUM_ATOMS; i+=gsz)
   {
      const int mid = material[i];
      const material_parameters_t mat = material_params[mid];

      const real_t alpha = mat.temperature_rescaling_alpha;
      const real_t sigma = mat.H_th_sigma;
      const real_t tc = mat.temperature_rescaling_Tc;
      const real_t norm_h = mat.applied_field_strength;

      const real_t resc_temp = (temp < tc) ? tc * POW(temp/tc, alpha) : temp;
      //const real_t resc_temp = select(temp,
      //                                tc*POW(temp/tc, alpha),
      //                                (uint_t)isless(temp, tc));
      const real_t sq_temp = sqrt(resc_temp);

      const size_t x = 3*i+0;
      const size_t y = 3*i+1;
      const size_t z = 3*i+2;

      const real_t3 grands = (real_t3)(gaussian_rand[x],
                                       gaussian_rand[y],
                                       gaussian_rand[z]);

      real_t3 field = sigma * sq_temp * grands;

      const real_t3 h = (real_t3)(mat.applied_field_unit_x,
                                  mat.applied_field_unit_y,
                                  mat.applied_field_unit_z);

      field += norm_h * h + Happ;

#ifdef USE_VECTOR_TYPE
      ext_field[i] = field + dip_field[i];
#else
      ext_field[x] = field.x + dip_field[x];
      ext_field[y] = field.y + dip_field[y];
      ext_field[z] = field.z + dip_field[z];
#endif
   }
}
