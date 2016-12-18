#include "cl_defs.h"

typedef struct
{
   real_t prefactor;
   real_t lambda_times_prefactor;
} heun_params_t;

__kernel
void llg_heun_predictor_step(const __global int *material_id,
                             __constant heun_params_t *heun_parameters,
                             __global real_t *x_spin,
                             __global real_t *y_spin,
                             __global real_t *z_spin,
                             const __global real_t *x_sp_field,
                             const __global real_t *y_sp_field,
                             const __global real_t *z_sp_field,
                             const __global real_t *x_ext_field,
                             const __global real_t *y_ext_field,
                             const __global real_t *z_ext_field,
                             __global real_t *dSx,
                             __global real_t *dSy,
                             __global real_t *dSz)
{
   const size_t gsz = get_global_size(0);

   for (size_t atom = get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      size_t mid = material_id[atom];

      real_t prefactor = heun_parameters[mid].prefactor;
      real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // initial spins
      real_t S_x = x_spin[atom];
      real_t S_y = y_spin[atom];
      real_t S_z = z_spin[atom];

      // total field
      real_t H_x = x_sp_field[atom] + x_ext_field[atom];
      real_t H_y = y_sp_field[atom] + y_ext_field[atom];
      real_t H_z = z_sp_field[atom] + z_ext_field[atom];

      // S cross H
      real_t SxH_x = S_y * H_z - S_z * H_y;
      real_t SxH_y = S_z * H_x - S_x * H_z;
      real_t SxH_z = S_x * H_y - S_y * H_x;

      // S cross (S cross H)
      real_t SxSxH_x = (S_y * SxH_z - S_z * SxH_y);
      real_t SxSxH_y = (S_z * SxH_x - S_x * SxH_z);
      real_t SxSxH_z = (S_x * SxH_y - S_y * SxH_x);

      // prefactor * (S cross H + lambda * S cross (S cross H))
      real_t DS_x = prefactor * SxH_x + lambdatpr * SxSxH_x;
      real_t DS_y = prefactor * SxH_y + lambdatpr * SxSxH_y;
      real_t DS_z = prefactor * SxH_z + lambdatpr * SxSxH_z;

      // change in S from predictor step
      dSx[atom] = DS_x;
      dSy[atom] = DS_y;
      dSz[atom] = DS_z;

      // predicted new spin
      real_t new_S_x = S_x + DS_x * DT;
      real_t new_S_y = S_y + DS_y * DT;
      real_t new_S_z = S_z + DS_z * DT;

      // normalization of spin
      real_t mod_s = RSQRT(
         new_S_x*new_S_x +
         new_S_y*new_S_y +
         new_S_z*new_S_z);

      x_spin[atom] = new_S_x * mod_s;
      y_spin[atom] = new_S_y * mod_s;
      z_spin[atom] = new_S_z * mod_s;
   }
}

__kernel
void llg_heun_corrector_step(const __global int *material_id,
                             __constant heun_params_t *hean_parameters,
                             __global real_t *x_spin,
                             __global real_t *y_spin,
                             __global real_t *z_spin,
                             const __global real_t *x_sp_field,
                             const __global real_t *y_sp_field,
                             const __global real_t *z_sp_field,
                             const __global real_t *x_ext_field,
                             const __global real_t *y_ext_field,
                             const __global real_t *z_ext_field,
                             const __global real_t *x_spin_buffer,
                             const __global real_t *y_spin_buffer,
                             const __global real_t *z_spin_buffer,
                             const __global real_t *dSx,
                             const __global real_t *dSy,
                             const __global real_t *dSz)
{
   size_t gsz = get_global_size(0);

   for (size_t atom = get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      size_t mid = material_id[atom];

      real_t prefactor = heun_parameters[mid].prefactor;
      real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // spins given by predictor step
      real_t Sx = x_spin[atom];
      real_t Sy = y_spin[atom];
      real_t Sz = z_spin[atom];

      // revised total field
      real_t H_x = x_sp_field[atom] + x_ext_field[atom];
      real_t H_y = y_sp_field[atom] + y_ext_field[atom];
      real_t H_z = z_sp_field[atom] + z_ext_field[atom];

      // S cross H
      real_t SxH_x = Sy * H_z - Sz * H_y;
      real_t SxH_y = Sz * H_x - Sx * H_z;
      real_t SxH_z = Sx * H_y - Sy * H_x;

      // S cross (S cross H)
      real_t SxSxH_x = Sy * SxH_z - Sz * SxH_y;
      real_t SxSxH_y = Sz * SxH_x - Sx * SxH_z;
      real_t SxSxH_z = Sx * SxH_y - Sy * SxH_x;

      // new change in S
      real_t DS_prime_x = prefactor * SxH_x + lambdatpr * SxSxH_x;
      real_t DS_prime_y = prefactor * SxH_y + lambdatpr * SxSxH_y;
      real_t DS_prime_z = prefactor * SxH_z + lambdatpr * SxSxH_z;

      // Heun step, using predictor and corrector spin changes
      real_t S_x = x_spin_buffer[atom] + 0.5 * (dSx[atom] + DS_prime_x) * DT;
      real_t S_y = y_spin_buffer[atom] + 0.5 * (dSy[atom] + DS_prime_y) * DT;
      real_t S_z = z_spin_buffer[atom] + 0.5 * (dSz[atom] + DS_prime_z) * DT;

      // normalization of spin
      real_t mod_s = RSQRT(
         S_x*S_x +
         S_y*S_y +
         S_z*S_z);

      x_spin[atom] = S_x * mod_s;
      y_spin[atom] = S_y * mod_s;
      z_spin[atom] = S_z * mod_s;
   }
}
