#include "cl_defs.h"

#ifdef USE_VECTOR_TYPE
typedef real_t3 T;
#else
typedef real_t T;
#endif

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
                             __global T *const restrict spin,
                             const __global T *const restrict sp_field,
                             const __global T *const restrict ext_field,
                             __global T *const restrict dS)
{
   const size_t gsz=get_global_size(0);

   for (size_t atom=get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
#ifdef USE_VECTOR_TYPE
      const real_t3 S = spin[atom];
      const real_t3 H = sp_field[atom] + ext_field[atom];
#else
      const size_t x = 3*atom+0;
      const size_t y = 3*atom+1;
      const size_t z = 3*atom+2;

      // initial spins
      const real_t3 S = (real_t3)(spin[x], spin[y], spin[z]);

      // total field
      const real_t3 H =
         (real_t3)(sp_field[x], sp_field[y], sp_field[z]) +
         (real_t3)(ext_field[x], ext_field[y], ext_field[z]);
#endif

      const size_t mid = material_id[atom];

      const real_t prefactor = heun_parameters[mid].prefactor;
      const real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // S cross H
      const real_t3 SxH = cross_product(S, H);

      // S cross (S cross H)
      const real_t3 SxSxH = cross_product(S, SxH);

      // prefactor * (S cross H + lambda * S cross (S cross H))
      const real_t3 Schange = prefactor * SxH + lambdatpr * SxSxH;

      // change in S from predictor step
#ifdef USE_VECTOR_TYPE
      dS[atom] = Schange;
#else
      dS[x] = Schange.x;
      dS[y] = Schange.y;
      dS[z] = Schange.z;
#endif

      // predicted new spin
      real_t3 new_S = S + Schange * DT;

      // normalization of spin
      const real_t rmod_s = RSQRT(new_S.x*new_S.x +
                                  new_S.y*new_S.y +
                                  new_S.z*new_S.z);

      new_S *= rmod_s;

#ifdef USE_VECTOR_TYPE
      spin[atom] = new_S;
#else
      spin[x] = new_S.x;
      spin[y] = new_S.y;
      spin[z] = new_S.z;
#endif
   }
}

__kernel
void llg_heun_corrector_step(const __global int *const restrict material_id,
                             __constant heun_params_t *const restrict heun_parameters,
                             __global T *const restrict spin,
                             const __global T *const restrict sp_field,
                             const __global T *const restrict ext_field,
                             const __global T *const restrict spin_buffer,
                             const __global T *const restrict dS)
{
   size_t gsz = get_global_size(0);

   for (size_t atom=get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
#ifdef USE_VECTOR_TYPE
      real_t3 S = spin[atom];
      real_t3 H = sp_field[atom] + ext_field[atom];
#else
      const size_t x = 3*atom+0;
      const size_t y = 3*atom+1;
      const size_t z = 3*atom+2;

      // spins given by predictor step
      real_t3 S = (real_t3)(spin[x], spin[y], spin[z]);

      // revised total field
      const real_t3 H =
         (real_t3)(sp_field[x], sp_field[y], sp_field[z]) +
         (real_t3)(ext_field[x], ext_field[y], ext_field[z]);
#endif

      const size_t mid = material_id[atom];

      const real_t prefactor = heun_parameters[mid].prefactor;
      const real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // S cross H
      const real_t3 SxH = cross_product(S, H);

      // S cross (S cross H)
      const real_t3 SxSxH = cross_product(S, SxH);

      // new change in S
      const real_t3 dS_prime = prefactor * SxH + lambdatpr * SxSxH;

#ifdef USE_VECTOR_TYPE
      S = spin_buffer[atom] + 0.5 * (dS[atom] + dS_prime) * DT;
#else
      // Heun step, using predictor and corrector spin changes
      S = (real_t3)(spin_buffer[x] + 0.5 * (dS[x] + dS_prime.x) * DT,
                    spin_buffer[y] + 0.5 * (dS[y] + dS_prime.y) * DT,
                    spin_buffer[z] + 0.5 * (dS[z] + dS_prime.z) * DT);
#endif

      // normalization of spin
      const real_t rmod_s = RSQRT(S.x*S.x +
                                  S.y*S.y +
                                  S.z*S.z);
      S *= rmod_s;

#ifdef USE_VECTOR_TYPE
      spin[atom] = S;
#else
      spin[x] = S.x;
      spin[y] = S.y;
      spin[z] = S.z;
#endif
   }
}
