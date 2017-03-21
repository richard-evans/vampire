//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

#include "cl_defs.h"

typedef struct
{
   real_t prefactor;
   real_t lambda_times_prefactor;
} heun_params_t;

__kernel
void llg_heun_predictor_step(const __global int *const restrict material_id,
                             const __constant heun_params_t *const restrict heun_parameters,
                                   __global real_t *const restrict spin,
                             const __global real_t *const restrict sp_field,
                             const __global real_t *const restrict ext_field,
                                   __global real_t *const restrict dS)
{
   const size_t gsz=get_global_size(0);

   for (size_t atom=get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      // initial spins
      const real_t3 S = vload3(atom, spin);

      // total field
      const real_t3 H = vload3(atom, sp_field) + vload3(atom, ext_field);

      const size_t mid = material_id[atom];

      const real_t prefactor = heun_parameters[mid].prefactor;
      const real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // S cross H
      const real_t3 SxH = cross(S, H);

      // S cross (S cross H)
      const real_t3 SxSxH = cross(S, SxH);

      // prefactor * (S cross H + lambda * S cross (S cross H))
      const real_t3 Schange = prefactor * SxH + lambdatpr * SxSxH;

      // change in S from predictor step
      vstore3(Schange, atom, dS);

      // predicted new spin
      real_t3 new_S = S + Schange * (real_t)DT;

      // normalization of spin
      new_S = normalize(new_S);

      vstore3(new_S, atom, spin);
   }
}

__kernel
void llg_heun_corrector_step(const __global int *const restrict material_id,
                             const __constant heun_params_t *const restrict heun_parameters,
                                   __global real_t *const restrict spin,
                             const __global real_t *const restrict sp_field,
                             const __global real_t *const restrict ext_field,
                             const __global real_t *const restrict spin_buffer,
                             const __global real_t *const restrict dS)
{
   size_t gsz = get_global_size(0);

   for (size_t atom=get_global_id(0); atom<NUM_ATOMS; atom+=gsz)
   {
      // spins given by predictor step
      real_t3 S = vload3(atom, spin);

      // revised total field
      const real_t3 H = vload3(atom, sp_field) + vload3(atom, ext_field);

      const size_t mid = material_id[atom];

      const real_t prefactor = heun_parameters[mid].prefactor;
      const real_t lambdatpr = heun_parameters[mid].lambda_times_prefactor;

      // S cross H
      const real_t3 SxH = cross(S, H);

      // S cross (S cross H)
      const real_t3 SxSxH = cross(S, SxH);

      // new change in S
      const real_t3 dS_prime = prefactor * SxH + lambdatpr * SxSxH;

      // Heun step, using predictor and corrector spin changes
      S = vload3(atom, spin_buffer) + 0.5 * (vload3(atom, dS) + dS_prime) * DT;

      // normalization of spin
      S = normalize(S);

      vstore3(S, atom, spin);
   }
}
