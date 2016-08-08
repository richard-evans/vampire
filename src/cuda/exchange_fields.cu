#include "exchange_fields.hpp"

#include "atoms.hpp"

#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"


#include <vector>
#include <thrust/device_vector.h>
#include <cusp/csr_matrix.h>

int calculate_exchange_fields(int, int);

namespace cu = vcuda::internal;

namespace vcuda
{

   namespace internal
   {

      namespace exchange
      {


         // Cusparse handle to use the csrmv routines
         cusparseHandle_t  cusparse_handle;

         // Cusparse descriptor for exchange matrices
         cusparseMatDescr_t J_descr;

         bool exchange_initialised = false;

         bool J_isot_initialised = false;
         bool J_vect_initialised = false;
         bool J_tens_initialised = false;


         int initialise_exchange()
         {

            // Initialise the cuSparse handle
            cusparse_call( cusparseCreate( &cusparse_handle) );

            // Initialise the matrix description
            cusparseCreateMatDescr(&J_descr);
            cusparseSetMatType(J_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(J_descr,CUSPARSE_INDEX_BASE_ZERO);


            check_device_memory(__FILE__,__LINE__);

            // Temporary vectors on host to copy exchange values from
            // They have to be declared outside of case statement
            std::vector<double> Jxx_vals_h;
            std::vector<double> Jyy_vals_h;
            std::vector<double> Jzz_vals_h;

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_xx_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_yy_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );

            cusp::csr_matrix < int, cu::cu_real_t, cusp::host_memory > J_zz_matrix_h (
                  ::atoms::num_atoms,
                  ::atoms::num_atoms,
                  ::atoms::neighbour_list_array.size()
                  );

            J_xx_matrix_h.row_offsets[0] = 0.0;
            J_yy_matrix_h.row_offsets[0] = 0.0;
            J_zz_matrix_h.row_offsets[0] = 0.0;
            for (int atom = 0; atom < ::atoms::num_atoms; atom++) {
               J_xx_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
               J_yy_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
               J_zz_matrix_h.row_offsets[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;
            }
            for (int i = 0; i < ::atoms::neighbour_list_array.size(); i++) {
               J_xx_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
               J_yy_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
               J_zz_matrix_h.column_indices[i] = ::atoms::neighbour_list_array[i];
            }


            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  for (int i = 0; i < ::atoms::neighbour_list_array.size(); i++) {
                     int iid = ::atoms::neighbour_interaction_type_array[i];
                     J_xx_matrix_h.values[i] = - ::atoms::i_exchange_list[iid].Jij;
                  }

                  cusp::convert(J_xx_matrix_h, J_xx_mat_d);

                  J_isot_initialised = true;
                  check_cuda_errors(__FILE__,__LINE__);

                  break;

               case 1: // Vector

                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  // Copy J values from vampire exchange list to values list
                  for( int i = 0; i < ::atoms::neighbour_list_array.size(); i++)
                  {
                     int iid = ::atoms::neighbour_interaction_type_array[i];
                     J_xx_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[0];
                     J_yy_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[1];
                     J_zz_matrix_h.values[i] = - ::atoms::v_exchange_list[iid].Jij[2];
                  }

                  cusp::convert(J_xx_matrix_h, J_xx_mat_d);
                  cusp::convert(J_yy_matrix_h, J_yy_mat_d);
                  cusp::convert(J_zz_matrix_h, J_zz_mat_d);

                  J_vect_initialised = true;
                  check_cuda_errors(__FILE__,__LINE__);

                  break;

               case 2: // Tensor
                  break;
            }

            exchange_initialised = true;

            std::cout << cu::atoms::limits[0] << "\t\t" << cu::atoms::limits[1] << "\t\t" << cu::atoms::limits[::atoms::num_atoms] << std::endl;


            std::cout << "Made matrix" << std::endl;
            check_device_memory(__FILE__,__LINE__);
            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         int finalise_exchange()
         {
            check_cuda_errors(__FILE__,__LINE__);
            // de-allocate all the data used in the exchange
            if( exchange_initialised)
            {
               cusparseDestroy( cusparse_handle);
               cusparseDestroyMatDescr( J_descr);
            }
            return EXIT_SUCCESS;
         }

         int calculate_exchange_fields()
         {

            check_cuda_errors(__FILE__,__LINE__);

            int * rowptrs_dptr = thrust::raw_pointer_cast( cu::atoms::limits.data() );
            int * colinds_dptr = thrust::raw_pointer_cast( cu::atoms::neighbours.data() );
            cu_real_t * Jxx_vals_dptr = thrust::raw_pointer_cast( Jxx_vals_d.data() );
            cu_real_t * Jyy_vals_dptr = thrust::raw_pointer_cast( Jyy_vals_d.data() );
            cu_real_t * Jzz_vals_dptr = thrust::raw_pointer_cast( Jzz_vals_d.data() );

            cu_real_t * Sx_dptr = thrust::raw_pointer_cast( cu::atoms::x_spin_array.data());
            cu_real_t * Sy_dptr = thrust::raw_pointer_cast( cu::atoms::y_spin_array.data());
            cu_real_t * Sz_dptr = thrust::raw_pointer_cast( cu::atoms::z_spin_array.data());
            cu_real_t * Hx_dptr = thrust::raw_pointer_cast( cu::x_total_spin_field_array.data());
            cu_real_t * Hy_dptr = thrust::raw_pointer_cast( cu::y_total_spin_field_array.data());
            cu_real_t * Hz_dptr = thrust::raw_pointer_cast( cu::z_total_spin_field_array.data());


            cusp::array1d_view <cu_real_array_t::iterator> x_spin_view (
                  cu::atoms::x_spin_array.begin(),
                  cu::atoms::x_spin_array.end()
                  );
            cusp::array1d_view <cu_real_array_t::iterator> y_spin_view (
                  cu::atoms::y_spin_array.begin(),
                  cu::atoms::y_spin_array.end()
                  );
            cusp::array1d_view <cu_real_array_t::iterator> z_spin_view (
                  cu::atoms::z_spin_array.begin(),
                  cu::atoms::z_spin_array.end()
                  );

            cusp::array1d_view <cu_real_array_t::iterator> x_total_spin_field_view (
                  cu::x_total_spin_field_array.begin(),
                  cu::x_total_spin_field_array.end()
                  );
            cusp::array1d_view <cu_real_array_t::iterator> y_total_spin_field_view (
                  cu::y_total_spin_field_array.begin(),
                  cu::y_total_spin_field_array.end()
                  );
            cusp::array1d_view <cu_real_array_t::iterator> z_total_spin_field_view (
                  cu::z_total_spin_field_array.begin(),
                  cu::z_total_spin_field_array.end()
                  );

            thrust::identity<cu_real_t>   identity;
            thrust::multiplies<cu_real_t> combine;
            thrust::plus<cu_real_t>       reduce;

            // cusparse csrmv calculates y = alpha * OP(A) * x + beta * y
            // where alpha and beta are scalar constants
            cu_real_t alpha = 1.0;
            cu_real_t beta = 1.0;

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

                  // FIXME This maybe boosted
                  // It should keep the old values stored in the spin field
                  // Since Jxx = Jyy = Jzz only the Jxx array is used
                  cusp::generalized_spgemm(
                        J_xx_mat_d,
                        x_spin_view,
                        x_total_spin_field_view,
                        identity, combine, reduce
                        );

                  cusp::generalized_spgemm(
                        J_xx_mat_d,
                        y_spin_view,
                        y_total_spin_field_view,
                        identity, combine, reduce
                        );

                  cusp::generalized_spgemm(
                        J_xx_mat_d,
                        z_spin_view,
                        z_total_spin_field_view,
                        identity, combine, reduce
                        );

                  break;

               case 1: // Vector exchange

                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

                  // Compute the Hx = Jxx * Sx exchange fields
                  cusparseTcsrmv(
                        cusparse_handle,
                         CUSPARSE_OPERATION_NON_TRANSPOSE,
                         ::atoms::num_atoms,
                         ::atoms::num_atoms,
                         ::atoms::total_num_neighbours,
                         &alpha,
                         J_descr,
                         Jxx_vals_dptr,
                         rowptrs_dptr,
                         colinds_dptr,
                         Sx_dptr,
                         &beta,
                         Hx_dptr);

                  // Compute the Hy = Jyy * Sy exchange fields
                  cusparseTcsrmv(
                        cusparse_handle,
                         CUSPARSE_OPERATION_NON_TRANSPOSE,
                         ::atoms::num_atoms,
                         ::atoms::num_atoms,
                         ::atoms::total_num_neighbours,
                         &alpha,
                         J_descr,
                         Jyy_vals_dptr,
                         rowptrs_dptr,
                         colinds_dptr,
                         Sy_dptr,
                         &beta,
                         Hy_dptr);

                  // Compute the Hz = Jzz * Sz exchange fields
                  cusparseTcsrmv(
                        cusparse_handle,
                         CUSPARSE_OPERATION_NON_TRANSPOSE,
                         ::atoms::num_atoms,
                         ::atoms::num_atoms,
                         ::atoms::total_num_neighbours,
                         &alpha,
                         J_descr,
                         Jzz_vals_dptr,
                         rowptrs_dptr,
                         colinds_dptr,
                         Sz_dptr,
                         &beta,
                         Hz_dptr);

                  break;
            }

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }
      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda

