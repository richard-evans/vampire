// Performs Box-Muller transform on
// uniform random numbers in range [0,1]
// to produce Gaussian distributed random numbers
// with mean 0 and standard deviation 1

#define TWO_PI 2.0*3.14159265358979323846

#ifdef OPENCL_DP
typedef double real_t
#else
typedef float  real_t
#endif

// native functions can be much faster but may have less accuracy
#ifdef OPENCL_USE_NATIVE_FUNCTIONS
#define SQRT(x) native_sqrt(x)
#define  LOG(x) native_log(x)
#define  COS(x) native_cos(x)
#else
#define SQRT(x) sqrt(x)
#define  LOG(x) log(x)
#define  COS(x) cos(x)
#endif


__kernel
void BoxMullerTransform(const __global uint *urands,
                        __global real_t *grands)
{
   const size_t gsz = get_global_size(0);

   for (unsigned id=get_global_id(0); id<N/2; id+=gsz)
   {
      real_t u1 = urands[2*id+0] / real_t(0xFFFFFFFF);
      real_t u2 = urands[2*id+1] / real_t(0xFFFFFFFF);

      real_t r = SQRT(-2*LOG(u1));
      real_t costheta = COS(TWO_PI * u2);

      grands[2*id+0] = r * costheta;
      grands[2*id+1] = r * SQRT(1 - costheta*costheta);
   }
}
