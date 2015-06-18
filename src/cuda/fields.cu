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
            double * x_sp_field, double * y_sp_field, double * z_sp_field

            )
      {
         size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
         size_t mid = material[tid];

         double field_x = 0.0;
         double field_y = 0.0;
         double field_z = 0.0;

         double sx = x_spin[tid];
         double sy = y_spin[tid];
         double sz = z_spin[tid];

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

         x_sp_field[tid] += field_x;
         y_sp_field[tid] += field_y;
         z_sp_field[tid] += field_z;

      }

   } /* internal */
#endif
} /* cuda */

