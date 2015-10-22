#include "cuda_utils.hpp"
#include "exchange_fields.hpp"
#include "data.hpp"
#include "internal.hpp"

#ifdef CUDA
namespace cu = vcuda::internal;
#endif

namespace vcuda
{
#ifdef CUDA
   namespace internal
   {
      __global__ void update_non_exchange_spin_fields (
            double * x_spin, double * y_spin, double * z_spin,
            int * material,
            cu::material_parameters_t * material_params,
            double * x_sp_field, double * y_sp_field, double * z_sp_field,
            int n_atoms
            )
      {
         for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {

            int mid = material[i];

            double field_x = 0.0;
            double field_y = 0.0;
            double field_z = 0.0;

            double sx = x_spin[i];
            double sy = y_spin[i];
            double sz = z_spin[i];

            /*
             * Scalar anisotropy
             */
            double ku = material_params[mid].ku;
            field_z -= 2.0 * ku * sz;

            double ex = material_params[mid].anisotropy_unit_x;
            double ey = material_params[mid].anisotropy_unit_y;
            double ez = material_params[mid].anisotropy_unit_z;

            double sdote = sx * ex + sy * ey + sz * ez;
            double sdote3 = sdote * sdote * sdote;
            double sdote5 = sdote3 * sdote * sdote;

            /*
             * Spherical harmonics
             */

            double scale = 0.6666666666666667;

            double k2 = material_params[mid].sh2;
            double k4 = material_params[mid].sh4;
            double k6 = material_params[mid].sh6;

            double ek2 = k2 * 3.0 * sdote;
            double ek4 = k4 * 0.125 * (140.0 * sdote3 - 60.0 *sdote);
            double ek6 = k6 * 0.0625 * (
                  1386.0 * sdote5 - 1260.0 * sdote3 + 210.0 * sdote);

            field_x += scale * ex * (ek2 + ek4 + ek6);
            field_y += scale * ey * (ek2 + ek4 + ek6);
            field_z += scale * ez * (ek2 + ek4 + ek6);

            /*
             * Lattice anisotropy
             */

            /*
             * TODO: add the temperature dependence
             */

            /*
             * TODO: communicate every timestep
             */

            double k_latt = 2.0 * material_params[mid].k_latt;
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

            x_sp_field[i] += field_x;
            y_sp_field[i] += field_y;
            z_sp_field[i] += field_z;
         }
      }

      __global__ void update_external_fields_kernel (
            int *  material,
            int * cell,
            vcuda::internal::material_parameters_t * material_params,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            double * x_ext_field, double * y_ext_field, double * z_ext_field,
            curandState * rand_state,
            double global_temperature,
            double Hx_app,
            double Hy_app,
            double Hz_app,
            int n_atoms
            )
      {

         /*
          * Thread and material identification
          */

         int tid = blockIdx.x * blockDim.x + threadIdx.x;

         for ( int i = tid;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {

            int mid = material[i];
            cu::material_parameters_t mat = material_params[mid];

            double field_x = 0.0;
            double field_y = 0.0;
            double field_z = 0.0;

            /*
             * TODO: HAMR fields
             */

            /*
             * Thermal fields
             */

            //double temp = mat.temperature;
            double temp = global_temperature;
            double alpha = mat.temperature_rescaling_alpha;
            double sigma = mat.H_th_sigma;
            double tc = mat.temperature_rescaling_Tc;
            double resc_temp = (temp < tc) ? tc * pow(temp / tc, alpha) : temp;
            double sq_temp = sqrt(resc_temp);

            field_x += sigma * sq_temp * curand_normal_double (
                  rand_state + tid);
            field_y += sigma * sq_temp * curand_normal_double (
                  rand_state + tid);
            field_z += sigma * sq_temp * curand_normal_double (
                  rand_state + tid);

            /*
             * Applied field
             */

            double norm_h = mat.applied_field_strength;
            double hx = mat.applied_field_unit_x;
            double hy = mat.applied_field_unit_y;
            double hz = mat.applied_field_unit_z;

            field_x += norm_h * hx;
            field_y += norm_h * hy;
            field_z += norm_h * hz;

            field_x += Hx_app;
            field_y += Hy_app;
            field_z += Hz_app;

            /*
             * TODO: FMR fields?
             */

            /*
             * Dipolar fields
             */

            field_x += x_dip_field[i];
            field_y += y_dip_field[i];
            field_z += z_dip_field[i];

            /*
             * Write back to main memory
             */

            x_ext_field[i] += field_x;
            y_ext_field[i] += field_y;
            z_ext_field[i] += field_z;

         }
      }

      __global__ void update_cell_magnetization (
            double * x_spin, double * y_spin, double * z_spin,
            int * material, int * cell,
            cu::material_parameters_t * material_params,
            double * x_mag, double * y_mag, double * z_mag,
            int n_atoms
            )
      {
         /*
          * TODO: This is an supremely naïve implementation
          *       the number of cells can be as big as the number of atoms
          *       so might as well leave it like this
          */

         for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {
            int mid = material[i];
            int cid = cell[i];
            double mu_s = material_params[mid].mu_s_si;
            cu::atomicAdd(&x_mag[cid], x_spin[i] * mu_s);
            cu::atomicAdd(&y_mag[cid], y_spin[i] * mu_s);
            cu::atomicAdd(&z_mag[cid], z_spin[i] * mu_s);
         }
      }

      __global__ void update_dipolar_fields (
            double * x_mag, double * y_mag, double * z_mag,
            double * x_coord, double * y_coord, double * z_coord,
            double * volume,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            int n_cells
            )
      {
         for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
               i < n_cells;
               i += blockDim.x * gridDim.x)
         {
            double mx = x_mag[i];
            double my = y_mag[i];
            double mz = z_mag[i];
            double cx = x_coord[i];
            double cy = y_coord[i];
            double cz = z_coord[i];
            /*
             * Inverse volume from the number of atoms in macro-cell
             */
            double vol_prefac = - 4.0 * M_PI / (3.0 * volume[i]);
            double prefactor = 1.0e+23; // 1e-7/1e30

            double field_x = vol_prefac * mx;
            double field_y = vol_prefac * my;
            double field_z = vol_prefac * mz;

            for (int j = 0; j < n_cells; j++)
            {
               if (i == j) continue;
               double omx = x_mag[i];
               double omy = y_mag[i];
               double omz = z_mag[i];

               double dx = x_coord[j] - cx;
               double dy = y_coord[j] - cy;
               double dz = z_coord[j] - cz;

               double drij = 1.0 / sqrtf (dx * dx + dy * dy + dz * dz);
               double drij3 = drij * drij * drij;

               double sdote = (
                     omx * dx * drij +
                     omy * dy * drij +
                     omz * dz * drij);

               field_x += (3.0 * sdote * dx * drij - omx) * drij3;
               field_y += (3.0 * sdote * dy * drij - omy) * drij3;
               field_z += (3.0 * sdote * dz * drij - omz) * drij3;
            }

            x_dip_field[i] = prefactor * field_x;
            y_dip_field[i] = prefactor * field_y;
            z_dip_field[i] = prefactor * field_z;
         }
      }

      __global__ void update_atomistic_dipolar_fields (
            double * x_cell_field, double * y_cell_field, double * z_cell_field,
            double * x_dip_field, double * y_dip_field, double * z_dip_field,
            int * cell, int n_atoms
            )
      {
         /*
          * TODO: Warning extremely naïve data access pattern
          */
         for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
               i < n_atoms;
               i += blockDim.x * gridDim.x)
         {
            int cid = cell[i];
            x_dip_field[i] = x_cell_field[cid];
            y_dip_field[i] = y_cell_field[cid];
            z_dip_field[i] = z_cell_field[cid];
         }
      }

      void update_spin_fields ()
      {
         /*
          * Fill the field vectors with zero
          */

         thrust::fill(
               cu::x_total_spin_field_array.begin(),
               cu::x_total_spin_field_array.end(),
               0.0);
         thrust::fill(
               cu::y_total_spin_field_array.begin(),
               cu::y_total_spin_field_array.end(),
               0.0);
         thrust::fill(
               cu::z_total_spin_field_array.begin(),
               cu::z_total_spin_field_array.end(),
               0.0);

         /*
          * Find the addresses in the device address space
          */

         double * d_x_spin = thrust::raw_pointer_cast(
               cu::atoms::x_spin_array.data());
         double * d_y_spin = thrust::raw_pointer_cast(
               cu::atoms::y_spin_array.data());
         double * d_z_spin = thrust::raw_pointer_cast(
               cu::atoms::z_spin_array.data());

         int * d_materials =
            thrust::raw_pointer_cast(cu::atoms::type_array.data());

         cu::material_parameters_t * d_material_params =
            thrust::raw_pointer_cast (cu::mp::materials.data());

         double * d_x_spin_field = thrust::raw_pointer_cast(
               cu::x_total_spin_field_array.data());
         double * d_y_spin_field = thrust::raw_pointer_cast(
               cu::y_total_spin_field_array.data());
         double * d_z_spin_field = thrust::raw_pointer_cast(
               cu::z_total_spin_field_array.data());

         cu::update_non_exchange_spin_fields <<< cu::grid_size, cu::block_size >>> (
               d_x_spin, d_y_spin, d_y_spin,
               d_materials, d_material_params,
               d_x_spin_field, d_y_spin_field, d_z_spin_field,
               ::atoms::num_atoms);

         check_cuda_errors (__FILE__, __LINE__);

         cu::exchange::calculate_exchange_fields ();

         check_cuda_errors (__FILE__, __LINE__);
      }

      void update_external_fields ()
      {

         thrust::fill(
               cu::x_total_external_field_array.begin(),
               cu::x_total_external_field_array.end(),
               0.0);
         thrust::fill(
               cu::y_total_external_field_array.begin(),
               cu::y_total_external_field_array.end(),
               0.0);
         thrust::fill(
               cu::z_total_external_field_array.begin(),
               cu::z_total_external_field_array.end(),
               0.0);

         /*
          * Find the addresses in the device address space
          */

         int * d_materials =
            thrust::raw_pointer_cast(cu::atoms::type_array.data());

         cu::material_parameters_t * d_material_params =
            thrust::raw_pointer_cast (cu::mp::materials.data());

         int * d_cells =
            thrust::raw_pointer_cast(cu::atoms::cell_array.data());

         double * d_x_dip_field = thrust::raw_pointer_cast(
               cu::x_dipolar_field_array.data());
         double * d_y_dip_field = thrust::raw_pointer_cast(
               cu::y_dipolar_field_array.data());
         double * d_z_dip_field = thrust::raw_pointer_cast(
               cu::z_dipolar_field_array.data());

         double * d_x_ext_field = thrust::raw_pointer_cast(
               cu::x_total_external_field_array.data());
         double * d_y_ext_field = thrust::raw_pointer_cast(
               cu::y_total_external_field_array.data());
         double * d_z_ext_field = thrust::raw_pointer_cast(
               cu::z_total_external_field_array.data());

         cu::update_external_fields_kernel <<< cu::grid_size, cu::block_size >>> (
               d_materials,
               d_cells,
               d_material_params,
               d_x_dip_field, d_y_dip_field, d_z_dip_field,
               d_x_ext_field, d_y_ext_field, d_z_ext_field,
               cu::d_rand_state,
               ::sim::temperature,
               ::sim::H_vec[0]*::sim::H_applied,
               ::sim::H_vec[1]*::sim::H_applied,
               ::sim::H_vec[2]*::sim::H_applied,
               ::atoms::num_atoms);

         check_cuda_errors (__FILE__, __LINE__);
      }

      void update_dipolar_fields ()
      {
         /*
          * Check if an update is required
          */

         if (::sim::time == ::demag::update_time) return;
         if (::sim::time % ::demag::update_rate) return;

         ::demag::update_time = ::sim::time;

         update_cell_magnetizations ();

         check_cuda_errors (__FILE__, __LINE__);

         /*
          * Figure out addresses in device memory space
          */

         double * d_x_mag = thrust::raw_pointer_cast(
               cu::cells::x_mag_array.data());
         double * d_y_mag = thrust::raw_pointer_cast(
               cu::cells::y_mag_array.data());
         double * d_z_mag = thrust::raw_pointer_cast(
               cu::cells::z_mag_array.data());

         double * d_x_coord = thrust::raw_pointer_cast(
               cu::cells::x_coord_array.data());
         double * d_y_coord = thrust::raw_pointer_cast(
               cu::cells::y_coord_array.data());
         double * d_z_coord = thrust::raw_pointer_cast(
               cu::cells::z_coord_array.data());

         double * d_volume = thrust::raw_pointer_cast(
               cu::cells::volume_array.data());

         double * d_x_cell_field = thrust::raw_pointer_cast(
               cu::cells::x_field_array.data());
         double * d_y_cell_field = thrust::raw_pointer_cast(
               cu::cells::y_field_array.data());
         double * d_z_cell_field = thrust::raw_pointer_cast(
               cu::cells::z_field_array.data());

         /*
          * Update cell dipolar fields
          */

         update_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
               d_x_mag, d_y_mag, d_z_mag,
               d_x_coord, d_y_coord, d_z_coord,
               d_volume,
               d_x_cell_field, d_y_cell_field, d_z_cell_field,
               ::cells::num_cells
               );

         check_cuda_errors (__FILE__, __LINE__);

         /*
          * Update atomistic dipolar fields
          */

         int * d_cells =
            thrust::raw_pointer_cast(cu::atoms::cell_array.data());

         double * d_x_atom_field = thrust::raw_pointer_cast(
               cu::x_dipolar_field_array.data());
         double * d_y_atom_field = thrust::raw_pointer_cast(
               cu::y_dipolar_field_array.data());
         double * d_z_atom_field = thrust::raw_pointer_cast(
               cu::z_dipolar_field_array.data());

         update_atomistic_dipolar_fields <<< cu::grid_size, cu::block_size >>> (
               d_x_cell_field, d_y_cell_field, d_z_cell_field,
               d_x_atom_field, d_y_atom_field, d_z_atom_field,
               d_cells,
               ::atoms::num_atoms
               );

         check_cuda_errors (__FILE__, __LINE__);
      }

      void update_cell_magnetizations ()
      {
         double * d_x_spin = thrust::raw_pointer_cast(
               cu::atoms::x_spin_array.data());
         double * d_y_spin = thrust::raw_pointer_cast(
               cu::atoms::y_spin_array.data());
         double * d_z_spin = thrust::raw_pointer_cast(
               cu::atoms::z_spin_array.data());

         int * d_materials =
            thrust::raw_pointer_cast(cu::atoms::type_array.data());

         int * d_cells =
            thrust::raw_pointer_cast(cu::atoms::cell_array.data());

         cu::material_parameters_t * d_material_params =
            thrust::raw_pointer_cast (cu::mp::materials.data());

         double * d_x_mag = thrust::raw_pointer_cast(
               cu::cells::x_mag_array.data());
         double * d_y_mag = thrust::raw_pointer_cast(
               cu::cells::y_mag_array.data());
         double * d_z_mag = thrust::raw_pointer_cast(
               cu::cells::z_mag_array.data());

         /*
          * Update cell magnetizations
          */

         thrust::fill(
               cu::cells::x_mag_array.begin(),
               cu::cells::x_mag_array.end(),
               0.0);
         thrust::fill(
               cu::cells::y_mag_array.begin(),
               cu::cells::y_mag_array.end(),
               0.0);
         thrust::fill(
               cu::cells::z_mag_array.begin(),
               cu::cells::z_mag_array.end(),
               0.0);

         update_cell_magnetization <<< cu::grid_size, cu::block_size >>> (
               d_x_spin, d_y_spin, d_z_spin,
               d_materials, d_cells,
               d_material_params,
               d_x_mag, d_y_mag, d_z_mag,
               ::atoms::num_atoms
               );

         check_cuda_errors (__FILE__, __LINE__);
      }


   } /* internal */
#endif
} /* cuda */

