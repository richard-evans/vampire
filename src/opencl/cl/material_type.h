#ifndef VOPENCL_MATERIAL_TYPE_H_
#define VOPENCL_MATERIAL_TYPE_H_

#include "cl_defs.h"

typedef struct
{
   real_t alpha;
   real_t gamma_rel;
   real_t mu_s_si;
   real_t i_mu_s_si;
   real_t k_latt;
   real_t sh2;
   real_t sh4;
   real_t sh6;
   real_t ku;
   real_t anisotropy_unit_x;
   real_t anisotropy_unit_y;
   real_t anisotropy_unit_z;
   real_t applied_field_strength;
   real_t applied_field_unit_x;
   real_t applied_field_unit_y;
   real_t applied_field_unit_z;
   real_t Kc1_SI;
   real_t temperature;
   real_t temperature_rescaling_alpha;
   real_t temperature_rescaling_Tc;
   real_t H_th_sigma;
} material_parameters_t;

#endif // VOPENCL_MATERIAL_TYPE_H_
