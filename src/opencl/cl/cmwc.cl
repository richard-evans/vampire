// Complementary multiply with carry generator

__kernel
void cmwc(__global uint *Q)
{   
   const size_t gsz = get_global_size(0);
   const size_t gid = get_global_id(0);

   uint c = Q[(gid+1)%N];

   for (size_t id=gid; id<N; id+=gsz)
   {
      const ulong t = 18782 * Q[id] + c;
      c = t >> 32;
      uint x = t + c;

      if (x < c)
      {
         ++x;
         ++c;
      }

      Q[id] = 0xFFFFFFFE - x;
   }
}
