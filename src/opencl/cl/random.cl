//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) S R H Morris 2017. All rights reserved.
//
//-----------------------------------------------------------------------------

// xorshift* uniform prng with Box Muller transform

#include "cl_defs.h"

ulong2 xorshift(ulong2 x)
{
   x ^= x >> 12;
   x ^= x << 25;
   x ^= x >> 27;

   return x;
}

__kernel
void gen_grands(__global ulong2  *const restrict state,
                __global real_t2 *const restrict grands)
{
   const size_t gid = get_global_id(0);
   const size_t gsz = get_global_size(0);

   for (size_t i=gid; i<(3*NUM_ATOMS)/2; i+=gsz)
   {
      ulong2 s = xorshift(state[i]);

      state[i] = s;

      // u.x, u.y are between 0 and 1
      const ulong c = 0x2545F4914F6CDD1Dul;
      s *= c;
      const real_t2 u =
#ifdef OPENCL_DP
         convert_double2(s)
#else
         convert_float2(s)
#endif
         /(real_t2)0xFFFFFFFFFFFFFFFFul;

      const real_t r = SQRT(-2*LOG(u.x));

      // TODO: look into sincos(), sinpi(), cospi()
      const real_t costheta = COS(2*PI*u.y);
      const real_t sintheta = SIN(2*PI*u.y);

      real_t2 g = r * (real_t2)(costheta, sintheta);

      grands[i] = g;
   }
}
