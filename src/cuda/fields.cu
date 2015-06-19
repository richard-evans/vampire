#include "internal.hpp"

namespace cuda
{
#ifdef CUDA
   namespace internal
   {

      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            size_t * material,
            vcuda::internal::material_parameters_t * material_params,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            size_t n_atoms
            )
      {
         size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
         for ( size_t i = tid;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {

            size_t mid = material[i];

            double field_x = 0.0;
            double field_y = 0.0;
            double field_z = 0.0;

            double sx = x_spin[i];
            double sy = y_spin[i];
            double sz = z_spin[i];

            /*
             * Scalar anisotropy
             */
            double ku = material_params[mid].ku;
            field_z -= 2.0 * ku * sz;

            /*
             * Second order uniaxial anisotropy
             */

            double ku2 = 4.0 * material_params[mid].ku2;

            double ex = material_params[mid].anisotropy_unit_x;
            double ey = material_params[mid].anisotropy_unit_y;
            double ez = material_params[mid].anisotropy_unit_z;

            double sdote = sx * ex + sy * ey + sz * ez;
            double sdote3 = sdote * sdote * sdote;
            field_x -= ku2 * ex * sdote3;
            field_y -= ku2 * ey * sdote3;
            field_z -= ku2 * ez * sdote3;

            /*
             * Sixth order o¿uniaxial anisotropy
             */

            double ku3 = 6.0 * material_params[mid].ku3;
            double sdote5 = sdote3 * sdote * sdote;
            field_x -= ku3 * ex * sdote5;
            field_y -= ku3 * ey * sdote5;
            field_z -= ku3 * ez * sdote5;

            /*
             * Spherical harmonics
             */

            double scale = 0.6666666666666667;

            double mu_s_si = material_params[mid].mu_s_SI;
            double k2 = material_params[mid].sh2 / mu_s_si;
            double k4 = material_params[mid].sh4 / mu_s_si;
            double k6 = material_params[mid].sh6 / mu_s_si;

            double ek2 = k2 * 3.0 * sdote;
            double ek4 = k4 * 0.125 * (140.0 * sdote3 - 60.0 *sdote);
            double ek6 = k6 * 0.0625 * (
                  1386.0 * sdote5 - 1260.0 * sdote3 + 210.0 * sdote);

            field_x += scale * ex * (ek2 + ek4 + ek6);
            field_y += scale * ey * (ek2 + ek4 + ek6);
            field_z += scale * ez * (ek2 + ek4 + ek6);

            /*
             * Lattice anisotropy
             * TODO: add the temperature dependence
             */

            double k_latt = 2.0 * material_params[mid].Klatt_SI / mu_s_si;
            field_x -= k_latt * ex * sdote;
            field_y -= k_latt * ey * sdote;
            field_z -= k_latt * ez * sdote;

            /*
             * TODO: Surface anisotropy?
             */

            /*
             * TODO: Lagrange multipliers?
             */

            /*
             * Write back to main memory
             */

            x_sp_field[i] += field_x;
            y_sp_field[i] += field_y;
            z_sp_field[i] += field_z;
         }
      }

      __global__ void update_external_fields (
            size_t * material, size_t * cell,
            vcuda::internal::material_parameters_t * material_params,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            curandState * rand_state, size_t n_atoms
            )
      {

         /*
          * Thread and material identification
          */

         size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

         for ( size_t i = blockIdx.x * blockDim.x + threadIdx.x;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {

            size_t mid = material[i];
            vcuda::internal::material_parameters_t mat = material_params[mid];

            double field_x = 0.0;
            double field_y = 0.0;
            double field_z = 0.0;

            /*
             * TODO: HAMR fields
             */

            /*
             * Thermal fields
             */

            double temp = mat.temperature;
            double alpha = mat.temperature_rescaling_alpha;
            double sigma = mat.H_th_sigma;
            double tc = mat.temperature_rescaling_Tc;
            double resc_temp = (temp < tc) ? tc * pow(temp / tc, alpha) : temp;
            double sq_temp = sqrt(resc_temp);

            field_x += sigma * sq_temp * curand_normal_double (rand_state + tid);
            field_y += sigma * sq_temp * curand_normal_double (rand_state + tid);
            field_z += sigma * sq_temp * curand_normal_double (rand_state + tid);

            /*
             * Applied field
             */

            double norm_h = mat.applied_field_strength;
            double hx = mat.applied_field_unit_x;
            double hy = mat.applied_field_unit_y;
            double hz = mat.applied_field_unit_z;

            field_x += norm_h * hx;
            field_y += norm_h * hy;
            field_z += norm_h * hz;

            /*
             * TODO: FMR fields?
             */

            /*
             * Dipolar fields
             */

            field_x += x_dip_field[i];
            field_y += y_dip_field[i];
            field_z += z_dip_field[i];

            /*
             * Write back to main memory
             */

            x_ext_field[i] += field_x;
            y_ext_field[i] += field_y;
            z_ext_field[i] += field_z;

         }
      }

   } /* internal */
#endif
} /* cuda */

