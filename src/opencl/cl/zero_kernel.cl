#include "cl_defs.h"

__kernel
void zero_buffers(real_t *const restrict b1,
                  real_t *const restrict b2,
                  real_t *const restrict b3)
{
   const size_t gid = get_global_id(0);

   b1[gid] = 0.0;
   b2[gid] = 0.0;
   b3[gid] = 0.0;
}
