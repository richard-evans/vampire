//------------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
// Reviewed: Andrea Meo 2022
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "typedefs.hpp"

// Conditional compilation of all cuda code
#ifdef CUDA

// namespace aliasing for brevity
namespace cu = vcuda::internal;

// vampire cuda namespace
namespace vcuda{

// module internal namespace
namespace internal{

//------------------------------------------------------------------------------
// Host function to calculate neel fields using gpu kernel
//------------------------------------------------------------------------------
void update_rotational_fields()
{
   cu::update_rotational_fields_kernel <<<cu::grid_size, cu::block_size>>> (cu::atoms::d_materials, cu::mp::d_material_params,
         cu::mp::d_material_params_r,
         cu::atoms::d_x_spin, cu::atoms::d_y_spin, cu::atoms::d_z_spin,
         cu::d_x_spin_field, cu::d_y_spin_field, cu::d_z_spin_field,
         ::atoms::num_atoms);

   check_cuda_errors(__FILE__, __LINE__);
}

__global__ void update_rotational_fields_kernel (
   int * material,
   cu::material_parameters_t * material_params,
   cu::material_parameters_rotational_t * material_params_r,
   cu_real_t * x_spin, cu_real_t * y_spin, cu_real_t * z_spin,
   cu_real_t * x_sp_field, cu_real_t * y_sp_field, cu_real_t * z_sp_field,
   int n_atoms
   )
{
   // TODO: Figure out why the results are still different from the serial output
   //
   for ( int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < n_atoms;
         i += blockDim.x * gridDim.x) {
      int mid = material[i];
      //Load parameters from memory
      cu::material_parameters_t material = material_params[mid];
      cu::material_parameters_rotational_t material_r = material_params_r[mid];

      // Initialise register to hold total spin field
      cu_real_t field_x = 0.0;
      cu_real_t field_y = 0.0;
      cu_real_t field_z = 0.0;

      // Registers to hold parts of spin field calculations
      cu_real_t x_component = 0.0;
      cu_real_t y_component = 0.0;
      cu_real_t z_component = 0.0;

      const cu_real_t sx = x_spin[i];
      const cu_real_t sy = y_spin[i];
      const cu_real_t sz = z_spin[i];

      const cu_real_t ex = material.anisotropy_unit_x;
      const cu_real_t ey = material.anisotropy_unit_y;
      const cu_real_t ez = material.anisotropy_unit_z;

      const cu_real_t fx = material_r.rotational_anisotropy_unit_x;
      const cu_real_t fy = material_r.rotational_anisotropy_unit_y;
      const cu_real_t fz = material_r.rotational_anisotropy_unit_z;

      const cu_real_t gx = material_r.last_axis_anisotropy_unit_x;
      const cu_real_t gy = material_r.last_axis_anisotropy_unit_y;
      const cu_real_t gz = material_r.last_axis_anisotropy_unit_z;

      const cu_real_t Sx = sx * fx + sy * fy + sz * fz;
      const cu_real_t Sy = sx * gx + sy * gy + sz * gz;
      const cu_real_t Sz = sx * ex + sy * ey + sz * ez;

      const cu_real_t Sx2 = Sx * Sx;
      const cu_real_t Sy2 = Sy * Sy;
      const cu_real_t Sz2 = Sz * Sz;

      const cu_real_t Sx3 = Sx2 * Sx;
      const cu_real_t Sy3 = Sy2 * Sy;

      const cu_real_t Sx4 = Sx2 * Sx2;
      const cu_real_t Sy4 = Sy2 * Sy2;
      const cu_real_t Sz4 = Sz2 * Sz2;

      const cu_real_t Sx2pSy2 = Sx2 + Sy2;

      //rotational_order_2_1
      const cu_real_t k2r1 = material_r.k2r1;

      field_x += (k2r1 * Sz * fx) + (k2r1 * Sx * ex);
      field_y += (k2r1 * Sz * fy) + (k2r1 * Sx * ey);
      field_z += (k2r1 * Sz * fz) + (k2r1 * Sx * ez);

      //rotational_order_2_1_odd
      const cu_real_t k2r1_odd = material_r.k2r1_odd;

      field_x += (k2r1_odd * Sz * gx) + (k2r1_odd * Sy * ex);
      field_y += (k2r1_odd * Sz * gy) + (k2r1_odd * Sy * ey);
      field_z += (k2r1_odd * Sz * gz) + (k2r1_odd * Sy * ez);

      //rotational_order_2_2
      const cu_real_t twok2r2 = 2.0 * material_r.k2r2;

      field_x += twok2r2 * Sx * fx;
      field_y += twok2r2 * Sx * fy;
      field_z += twok2r2 * Sx * fz;

      field_x -= twok2r2 * Sy * gx;
      field_y -= twok2r2 * Sy * gy;
      field_z -= twok2r2 * Sy * gz;

      //rotational_order_2_2_odd
      const cu_real_t twok2r2_odd = 2.0 * material_r.k2r2_odd;

      field_x += (twok2r2_odd * Sy * fx) + (twok2r2_odd * Sx * gx);
      field_y += (twok2r2_odd * Sy * fy) + (twok2r2_odd * Sx * gy);
      field_z += (twok2r2_odd * Sy * fz) + (twok2r2_odd * Sx * gz);

      //rotational_order_4_1
      const cu_real_t k4r1 = material_r.k4r1;

      x_component = k4r1 * Sz * (Sz2 - (3.0/7.0));
      z_component = k4r1 * Sx * (3.0 * Sz2 - (3.0/7.0));

      field_x += x_component * fx + z_component * ex;
      field_y += x_component * fy + z_component * ey;
      field_z += x_component * fz + z_component * ez;

      //rotational_order_4_1_odd
      const cu_real_t k4r1_odd = material_r.k4r1_odd;

      y_component = k4r1_odd * Sz * (Sz2 - (3.0/7.0));
      z_component = k4r1_odd * Sy * (3.0 * Sz2 - (3.0/7.0));

      field_x += y_component * gx + z_component * ex;
      field_y += y_component * gy + z_component * ey;
      field_z += y_component * gz + z_component * ez;

      //rotational_order_4_2
      const cu_real_t k4r2 = material_r.k4r2;

      x_component = k4r2 * ((Sx3 * 4.0) - (Sx * (12.0/7.0)));
      y_component = k4r2 * ((Sy3 * 4.0) - (Sy * (12.0/7.0)));

      field_x += - x_component * fx + y_component * gx;
      field_y += - x_component * fy + y_component * gy;
      field_z += - x_component * fz + y_component * gz;

      //rotational_order_4_2_odd
      const cu_real_t k4r2_odd = material_r.k4r2_odd;

      x_component = 2.0 * k4r2_odd * Sy * ((6.0/7.0) - 3.0 * Sx2 - Sy2);
      y_component = 2.0 * k4r2_odd * Sx * ((6.0/7.0) - 3.0 * Sy2 - Sx2);

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_4_3
      const cu_real_t k4r3 = material_r.k4r3;

      x_component = k4r3 * 3.0 * Sz * (Sx2 - Sy2);
      y_component = - k4r3 * 6.0 * Sx * Sy * Sz;
      z_component = k4r3 * Sx * (Sx2 - 3.0 * Sy2);

      field_x += x_component * fx + y_component * gx + z_component * ex;
      field_y += x_component * fy + y_component * gy + z_component * ey;
      field_z += x_component * fz + y_component * gz + z_component * ez;

      //rotational_order_4_3_odd
      const cu_real_t k4r3_odd = material_r.k4r3_odd;

      x_component = k4r3_odd * 6.0 * Sx * Sy * Sz;
      y_component = k4r3_odd * 3.0 * Sz * (Sx2 - Sy2);
      z_component = k4r3_odd * Sy * (3.0 * Sx2 - Sy2);

      field_x += x_component * fx + y_component * gx + z_component * ex;
      field_y += x_component * fy + y_component * gy + z_component * ey;
      field_z += x_component * fz + y_component * gz + z_component * ez;

      //rotational_order_4_4
      const cu_real_t k4r4 = material_r.k4r4;

      x_component = 4.0 * k4r4 * Sx * (Sx2 - 3.0 * Sy2);
      y_component = 4.0 * k4r4 * Sy * (Sy2 - 3.0 * Sx2);

      field_x += x_component * fx;
      field_y += x_component * fy;
      field_z += x_component * fz;

      field_x += y_component * gx;
      field_y += y_component * gy;
      field_z += y_component * gz;

      //rotational_order_4_4_odd
      const cu_real_t k4r4_odd = material_r.k4r4_odd;

      x_component = 4.0 * k4r4_odd * Sy * (3.0 * Sx2 - Sy2);
      y_component = 4.0 * k4r4_odd * Sx * (Sx2 - 3.0 * Sy2);

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_6_1
      const cu_real_t k6r1 = material_r.k6r1;

      x_component = k6r1 * Sz * (Sz4 - (10.0/11.0) * Sz2 + (5.0/33.0));
      z_component = k6r1 * Sx * 5.0 * (Sz4 - (6.0/11.0) * Sz2 + (1.0/33.0));

      field_x += z_component * ex + x_component * fx;
      field_y += z_component * ey + x_component * fy;
      field_z += z_component * ez + x_component * fz;
      //rotational_order_6_1_odd
      const cu_real_t k6r1_odd = material_r.k6r1_odd;

      y_component = k6r1_odd * Sz * (Sz4 - (10.0/11.0) * Sz2 + (5.0/33.0));
      z_component = k6r1_odd * Sy * 5.0 * (Sz4 - (6.0/11.0) * Sz2 + (1.0/33.0));

      field_x += z_component * ex + y_component * gx;
      field_y += z_component * ey + y_component * gy;
      field_z += z_component * ez + y_component * gz;

      //rotational_order_6_2
      const cu_real_t two_k6r2 = 2.0 * material_r.k6r2;

      x_component = two_k6r2 * Sx * (Sx2pSy2 * (3.0 * Sx2 - Sy2) + (16.0/33.0) * (1.0 - 6.0 * Sx2));
      y_component = two_k6r2 * Sy * (Sx2pSy2 * (Sx2 - 3.0 * Sy2) + (16.0/33.0) * (6.0 * Sy2 - 1.0));

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_6_2_odd
      const cu_real_t two_k6r2_odd = 2.0 * material_r.k6r2_odd;

      x_component = two_k6r2_odd * Sy * (5.0 * Sx4 + 6.0 * (Sx2 * Sy2) + Sy4 - (16.0/11.0) * (3.0 * Sx2 + Sy2) + (16.0/33.0));
      y_component = two_k6r2_odd * Sx * (Sx4 + 6.0 * (Sx2 * Sy2) + 5.0 * Sy4 - (16.0/11.0) * (Sx2 + 3.0 * Sy2) + (16.0/33.0));

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_6_3
      const cu_real_t k6r3 = material_r.k6r3;

      x_component = k6r3 * Sz * 3.0 * (Sz4 + 4.0 * (Sx2 * Sz2) - (1.0/11.0) * (14.0 * Sz2 + 12.0 * Sx2 - 3.0));
      z_component = k6r3 * Sx * 3.0 * (5.0 * Sz4 + 4.0 * (Sx2 * Sz2) - (1.0/11.0) * (42.0 * Sz2 + 4.0 * Sx2 - 3.0));

      field_x += z_component * ex + x_component * fx;
      field_y += z_component * ey + x_component * fy;
      field_z += z_component * ez + x_component * fz;

      //rotational_order_6_3_odd
      const cu_real_t k6r3_odd = material_r.k6r3_odd;

      y_component = - k6r3_odd * Sz * 3.0 * (Sz4 + 4.0 * (Sy2*Sz2) - (1.0/11.0) * (14.0 * Sz2 + 12.0 * Sy2 - 3.0));
      z_component = - k6r3_odd * Sy * 3.0 * (5.0 * Sz4 + 4.0 * (Sy2*Sz2) - (1.0/11.0) * (42.0 * Sz2 + 4.0 * Sy2 - 3.0));

      field_x += z_component * ex + y_component * gx;
      field_y += z_component * ey + y_component * gy;
      field_z += z_component * ez + y_component * gz;

      //rotational_order_6_4
      const cu_real_t two_k6r4 = 2.0 * material_r.k6r4;

      x_component = two_k6r4 * Sx * (3.0 * Sx4 - 10.0 * (Sx2 * Sy2) - 5.0 * Sy4 - (20.0/11.0) * (Sx2 - 3.0 * Sy2));
      y_component = two_k6r4 * Sy * (3.0 * Sy4 - 10.0 * (Sx2 * Sy2) - 5.0 * Sx4 - (20.0/11.0) * (Sy2 - 3.0 * Sx2));

      field_x -= x_component * fx + y_component * gx;
      field_y -= x_component * fy + y_component * gy;
      field_z -= x_component * fz + y_component * gz;

      //rotational_order_6_4_odd
      const cu_real_t four_k6r4_odd = 4.0 * material_r.k6r4_odd;

      x_component = four_k6r4_odd * Sy * (Sy4 - 5.0 * Sx4 + (10.0/11.0) * (3.0 * Sx2 - Sy2));
      y_component = four_k6r4_odd * Sx * (5.0 * Sy4 - Sx4 + (10.0/11.0) * (Sx2 - 3.0 * Sy2));

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_6_5
      const cu_real_t k6r5 = material_r.k6r5;

      x_component = k6r5 * Sz * 5.0 * (16.0 * Sx4 + 12.0 * (Sx2 * Sz2) + Sz4 - 2.0 * (6.0 * Sx2 + Sz2) + 1.0);
      z_component = k6r5 * Sx * (16.0 * Sx4 + 60.0 * (Sx2 * Sz2) + 25.0 * Sz4 - 10.0 * (2.0 * Sx2 + 3.0 * Sz2) + 5.0);

      field_x += z_component * ex + x_component * fx;
      field_y += z_component * ey + x_component * fy;
      field_z += z_component * ez + x_component * fz;

      //rotational_order_6_5_odd
      const cu_real_t k6r5_odd = material_r.k6r5_odd;

      y_component = k6r5_odd * Sz * 5.0 * (16.0 * Sy4 + 12.0 * (Sy2 * Sz2) + Sz4 - 2.0 * (6.0 * Sy2 + Sz2) + 1.0);
      z_component = k6r5_odd * Sy * (16.0 * Sy4 + 60.0 * (Sy2 * Sz2) + 25.0 * Sz4 - 10.0 * (2.0 * Sy2 + 3.0 * Sz2) + 5.0);

      field_x += z_component * ex + y_component * gx;
      field_y += z_component * ey + y_component * gy;
      field_z += z_component * ez + y_component * gz;

      //rotational_order_6_6
      const cu_real_t six_k6r6 = 6.0 * material_r.k6r6;

      x_component = six_k6r6 * Sx * (Sx4 - 10.0 * (Sx2 * Sy2) + 5.0 * Sy4);
      y_component = - six_k6r6 * Sy * (Sy4 - 10.0 * (Sx2 * Sy2) + 5.0 * Sx4);

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //rotational_order_6_6_odd
      const cu_real_t six_k6r6_odd = 6.0 * material_r.k6r6_odd;

      x_component = six_k6r6_odd * Sy * (5.0 * Sx4 - 10.0 * (Sx2 * Sy2) + Sy4);
      y_component = six_k6r6_odd * Sx * (Sx4 - 10.0 * (Sx2 * Sy2) + 5.0 * Sy4);

      field_x += x_component * fx + y_component * gx;
      field_y += x_component * fy + y_component * gy;
      field_z += x_component * fz + y_component * gz;

      //Apply all effects to the spin field
      x_sp_field[i] += field_x;
      y_sp_field[i] += field_y;
      z_sp_field[i] += field_z;
   }
}

} // end of internal namespace

} // end of vcuda namespace

#endif
