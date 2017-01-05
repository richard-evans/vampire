#include "cl_defs.h"

typedef struct
{
   real_t prefactor;
   real_t lambda_times_prefactor;
} heun_params_t;

real_t3 cross_product(const real_t3 A, const real_t3 B)
{
   return (real_t3)(A.y*B.z - A.z*B.y,
                    A.z*B.x - A.x*B.z,
                    A.x*B.y - A.y*B.x);
}

__kernel
void llg_heun_predictor_step(const __global int *const restrict material_id,
                             __constant heun_params_t *const restrict heun_parameters,
                             __global real_t3 *const restrict spin,
                             const __global real_t3 *const restrict sp_field,
                             const __global real_t3 *const restrict ext_field,
                             __global real_t3 *const restrict dS)
{
   const size_t gsz = get_global_size(0);

   for (size_t atom = get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      size_t mid = material_id[atom];

      real_t prefactor = heun_parameters[mid].prefactor;
      real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // initial spins
      real_t3 S = spin[atom];

      // total field
      real_t3 H = sp_field[atom] + ext_field[atom];

      // S cross H
      real_t3 SxH = cross_product(S, H);

      // S cross (S cross H)
      real_t3 SxSxH = cross_product(S, SxH);

      // prefactor * (S cross H + lambda * S cross (S cross H))
      real_t3 DS = prefactor * SxH + lambdatpr * SxSxH;

      // change in S from predictor step
      dS[atom] = DS;

      // predicted new spin
      real_t3 new_S = S + DS * DT;

      // normalization of spin
      real_t mod_s = RSQRT(new_S.x*new_S.x +
                           new_S.y*new_S.y +
                           new_S.z*new_S.z);

      spin[atom] = new_S * mod_s;
   }
}

__kernel
void llg_heun_corrector_step(const __global int *const restrict material_id,
                             __constant heun_params_t *const restrict heun_parameters,
                             __global real_t3 *const restrict spin,
                             const __global real_t3 *const restrict sp_field,
                             const __global real_t3 *const restrict ext_field,
                             const __global real_t3 *const restrict spin_buffer,
                             const __global real_t3 *const restrict dS)
{
   size_t gsz = get_global_size(0);

   for (size_t atom = get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      size_t mid = material_id[atom];

      real_t prefactor = heun_parameters[mid].prefactor;
      real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // spins given by predictor step
      real_t3 S = spin[atom];

      // revised total field
      real_t3 H = sp_field + ext_field;

      // S cross H
      real_t3 SxH = cross_product(S, H);

      // S cross (S cross H)
      real_t3 SxSxH = cross_product(S, SxH);

      // new change in S
      real_t3 dS_prime = prefactor * SxH + lambdatpr * SxSxH;

      // Heun step, using predictor and corrector spin changes
      S = spin_buffer[atom] + 0.5 * (dS[atom] + dS_prime) * DT;

      // normalization of spin
      real_t mod_s = RSQRT(S.x*S.x +
                           S.y*S.y +
                           S.z*S.z);

      spin[atom] = S * mod_s;
   }
}
