#include "exchange_fields.hpp"

#include "atoms.hpp"

#include "cuda_utils.hpp"
#include "internal.hpp"
#include "data.hpp"


#include <vector>
#include <thrust/device_vector.h>

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

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  // Copy J values from vampire exchange list to values list
                  for( int i = 0; i < ::atoms::neighbour_list_array.size(); i++)
                  {
                     int iid = ::atoms::neighbour_interaction_type_array[i]; // interaction id
                     double Jij= ::atoms::i_exchange_list[iid].Jij;

                     // -ve required for convention
                     Jxx_vals_h.push_back(-Jij);

                  }

                  Jxx_vals_d.resize( Jxx_vals_h.size() );

                  thrust::copy(
                        Jxx_vals_h.begin(),
                        Jxx_vals_h.end(),
                        Jxx_vals_d.begin()
                        );

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
                     int iid = ::atoms::neighbour_interaction_type_array[i]; // interaction id
                     Jxx_vals_h.push_back( -::atoms::v_exchange_list[iid].Jij[0]);
                     Jyy_vals_h.push_back( -::atoms::v_exchange_list[iid].Jij[1]);
                     Jzz_vals_h.push_back( -::atoms::v_exchange_list[iid].Jij[2]);

                  }

                  thrust::copy(
                        Jxx_vals_h.begin(),
                        Jxx_vals_h.end(),
                        Jxx_vals_d.begin()
                        );
                  thrust::copy(
                        Jyy_vals_h.begin(),
                        Jyy_vals_h.end(),
                        Jyy_vals_d.begin()
                        );
                  thrust::copy(
                        Jyy_vals_h.begin(),
                        Jyy_vals_h.end(),
                        Jyy_vals_d.begin()
                        );

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
            double * Jxx_vals_dptr = thrust::raw_pointer_cast( Jxx_vals_d.data() );
            double * Jyy_vals_dptr = thrust::raw_pointer_cast( Jyy_vals_d.data() );
            double * Jzz_vals_dptr = thrust::raw_pointer_cast( Jzz_vals_d.data() );

            double * Sx_dptr = thrust::raw_pointer_cast( cu::atoms::x_spin_array.data());
            double * Sy_dptr = thrust::raw_pointer_cast( cu::atoms::y_spin_array.data());
            double * Sz_dptr = thrust::raw_pointer_cast( cu::atoms::z_spin_array.data());
            double * Hx_dptr = thrust::raw_pointer_cast( cu::x_total_spin_field_array.data());
            double * Hy_dptr = thrust::raw_pointer_cast( cu::y_total_spin_field_array.data());
            double * Hz_dptr = thrust::raw_pointer_cast( cu::z_total_spin_field_array.data());


            // cusparse csrmv calculates y = alpha * OP(A) * x + beta * y
            // where alpha and beta are scalar constants
            double alpha = 1.0;
            double beta = 1.0;

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();


                  // Since Jxx = Jyy = Jzz only the Jxx array is used
                  cusparseDcsrmv(
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

                  cusparseDcsrmv(
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
                         Sy_dptr,
                         &beta,
                         Hy_dptr);

                  cusparseDcsrmv(
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
                         Sz_dptr,
                         &beta,
                         Hz_dptr);

                  break;

               case 1: // Vector exchange

                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

                  // Compute the Hx = Jxx * Sx exchange fields
                  cusparseDcsrmv(
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
                  cusparseDcsrmv(
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
                  cusparseDcsrmv(
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

