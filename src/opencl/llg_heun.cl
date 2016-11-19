#ifdef OPENCL_DP
typedef double real_t
#else
typedef float  real_t
#endif

struct heun_params_t
{
   real_t prefactor;
   real_t lambda_times_prefactor;
};

__kernel
void llg_heun_predictor_step(const __global int *material_id,
                             const __global heun_params_t *hean_params,
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
                             __global real_t *dSz,
                             const real_t dt,
                             const size_t num_atoms)
{
   size_t gsz = get_global_size(0);

   for (size_t atom = get_global_id(0); atom<num_atoms; atoms+=gsz)
   {
      size_t mid = material_id[atom];

      real_t prefactor = heun_parameters[mid].prefactor;
      real_t lambdatpr = heun_parameters[mid].lambda_time_prefactor;

      real_t sx = x_spin[atom];
      real_t sy = y_spin[atom];
      real_t sz = z_spin[atom];

      real_t H_x = x_sp_field[atom] + x_ext_field[atom];
      real_t H_y = y_sp_field[atom] + y_ext_field[atom];
      real_t H_z = z_sp_field[atom] + z_ext_field[atom];

      // s cross h
      real_t sxh_x = sy * H_z - sz * H_y;
      real_t sxh_y = sz * H_x - sx * H_z;
      real_t sxh_z = sx * H_y - sy * H_x;

      real_t Ds_x = prefactor * sxh_x + lambdatpr * (sy * sxh_z - sz * sxh_y);
      real_t Ds_y = prefactor * sxh_y + lambdatpr * (sz * sxh_x - sx * sxh_z);
      real_t Ds_z = prefactor * sxh_z + lambdatpr * (sx * sxh_y - sy * sxh_x);

      dSx[atom] = Ds_x;
      dSy[atom] = Ds_y;
      dSz[atom] = Ds_z;

      real_t new_spin_x = sx + Ds_x * dt;
      real_t new_spin_y = sy + Ds_y * dt;
      real_t new_spin_z = sz + Ds_z * dt;

      real_t mod_s = rsqrt(
         new_spin_x*new_spin_x +
         new_spin_y*new_spin_y +
         new_spin_z*new_spin_z);

      x_spin[atom] = new_spin_x * mod_s;
      y_spin[atom] = new_spin_y * mod_s;
      z_spin[atom] = new_spin_z * mod_s;
   }
}
