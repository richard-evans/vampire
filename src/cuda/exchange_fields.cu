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

         bool exchange_initialised = false;

         cusparseHandle_t  cusparse_handle;

         // Exchange matrix stored as CSR sparse matrix
         thrust::device_vector<int> rowptrs_d;
         thrust::device_vector<int> colinds_d;

         thrust::device_vector<double> Jxx_vals_d;
         thrust::device_vector<double> Jyy_vals_d;
         thrust::device_vector<double> Jzz_vals_d;

         bool Jxx_initialised = false;

         cusparseMatDescr_t J_descr;

         int initialise_exchange()
         {

            // Initialise the cuSparse handle
            cusparse_call( cusparseCreate( &cusparse_handle) );

            cusparseCreateMatDescr(&J_descr);
            cusparseSetMatType(J_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(J_descr,CUSPARSE_INDEX_BASE_ZERO);

            size_t avail;
            size_t total;
            cudaMemGetInfo( &avail, &total );
            size_t used = total - avail;
            std::cout << "Total: " << total << " , Used: " << used << " , Free: " << avail << std::endl;


            // Transfer the row ptrs and col indices to the device
            std::vector<int> rowptrs_h( ::atoms::num_atoms + 1, 0);
            for( int atom = 0; atom < ::atoms::num_atoms; atom++)
               rowptrs_h[atom+1] = ::atoms::neighbour_list_end_index[atom]+1;

            rowptrs_d.resize( ::atoms::num_atoms + 1);
            colinds_d.resize( ::atoms::neighbour_list_array.size() );


            thrust::copy(
                  rowptrs_h.begin(),
                  rowptrs_h.end(),
                  rowptrs_d.begin()
                  );

            thrust::copy(
                  ::atoms::neighbour_list_array.begin(),
                  ::atoms::neighbour_list_array.end(),
                  colinds_d.begin()
                  );

            std::vector<double> Jxx_vals_h( ::atoms::neighbour_list_array.size());
            std::vector<double> Jyy_vals_h( ::atoms::neighbour_list_array.size());
            std::vector<double> Jzz_vals_h( ::atoms::neighbour_list_array.size());

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  std::cerr << "Setting up isotropic exchange interaction" << std::endl;

                  for( int i = 0; i < ::atoms::neighbour_list_array.size(); i++)
                  {
                     int iid = ::atoms::neighbour_interaction_type_array[i]; // interaction id
                     double Jij= ::atoms::i_exchange_list[iid].Jij;

                     Jxx_vals_h[i] = -Jij;

                  }

                  Jxx_vals_d.resize( Jxx_vals_h.size() );

                  thrust::copy(
                        Jxx_vals_h.begin(),
                        Jxx_vals_h.end(),
                        Jxx_vals_d.begin()
                        );

                  Jxx_initialised = true;
                  check_cuda_errors(__FILE__,__LINE__);

                  break;

               case 1: // Vector

                  for( int i = 0; i < ::atoms::neighbour_list_array.size(); i++)
                  {
                     int iid = ::atoms::neighbour_interaction_type_array[i]; // interaction id
                     Jxx_vals_h[i] = -::atoms::v_exchange_list[iid].Jij[0];
                     Jyy_vals_h[i] = -::atoms::v_exchange_list[iid].Jij[1];
                     Jzz_vals_h[i] = -::atoms::v_exchange_list[iid].Jij[2];

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

                  Jxx_initialised = true;
                  check_cuda_errors(__FILE__,__LINE__);

                  break;

               case 2: // Tensor
                  break;
            }

            exchange_initialised = true;

            std::cout << "Made matrix" << std::endl;
            size_t avail_new, total_new;
            cudaError_t err = cudaMemGetInfo( &avail_new, &total_new );
            if( err != cudaSuccess)
            {
               std::cout << "Error fetching memory __LINE__" << std::endl;
            }
            size_t used_new = total_new - avail_new;
            std::cout << "Total: " << total_new << " , Used: " << used_new << " , Free: " << avail_new << std::endl;
            std::cout << "Memory used by matrix on device: " << used_new - used << " or in mb: " << float(used_new - used)/(1024.0*1024.0) << std::endl;

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }


         int finalise_exchange()
         {
            check_cuda_errors(__FILE__,__LINE__);
            if( exchange_initialised)
            {

               cusparseDestroy( cusparse_handle);
               cusparseDestroyMatDescr( J_descr);

            }

            // Manually call the destructors for the thrust vectors
            rowptrs_d.~device_vector();
            colinds_d.~device_vector();

            Jxx_vals_d.~device_vector();
            Jyy_vals_d.~device_vector();
            Jzz_vals_d.~device_vector();


            return EXIT_SUCCESS;
         }

         int calculate_exchange_fields()
         {

            check_cuda_errors(__FILE__,__LINE__);
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


            int * rowptrs_dptr = thrust::raw_pointer_cast( rowptrs_d.data() );
            int * colinds_dptr = thrust::raw_pointer_cast( colinds_d.data() );
            double * Jxx_vals_dptr = thrust::raw_pointer_cast( Jxx_vals_d.data() );
            double * Jyy_vals_dptr = thrust::raw_pointer_cast( Jyy_vals_d.data() );
            double * Jzz_vals_dptr = thrust::raw_pointer_cast( Jzz_vals_d.data() );

            double * Sx_dptr = thrust::raw_pointer_cast( cu::atoms::x_spin_array.data());
            double * Sy_dptr = thrust::raw_pointer_cast( cu::atoms::y_spin_array.data());
            double * Sz_dptr = thrust::raw_pointer_cast( cu::atoms::z_spin_array.data());
            double * Hx_dptr = thrust::raw_pointer_cast( cu::x_total_spin_field_array.data());
            double * Hy_dptr = thrust::raw_pointer_cast( cu::y_total_spin_field_array.data());
            double * Hz_dptr = thrust::raw_pointer_cast( cu::z_total_spin_field_array.data());

            double alpha = 1.0;
            double beta = 0.0;

            switch( ::atoms::exchange_type)
            {
               case 0: // Isotropic

                  //--------------------------------------------------------------
                  // Exchange is isotropic so Jxx = Jyy = Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

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

                  std::cout << "Vector Exchange" << std::endl;
                  //--------------------------------------------------------------
                  // Exchange is diagonal so Jxx != Jyy != Jzz
                  // and Jxy = Jxz = Jyx = 0
                  //--------------------------------------------------------------

                  if( !exchange_initialised) initialise_exchange();

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
                         Jyy_vals_dptr,
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
                         Jzz_vals_dptr,
                         rowptrs_dptr,
                         colinds_dptr,
                         Sz_dptr,
                         &beta,
                         Hz_dptr);

                  break;

            }

            ::atoms::x_total_spin_field_array[0] = 0.0;
            ::atoms::y_total_spin_field_array[0] = 0.0;
            ::atoms::z_total_spin_field_array[0] = 0.0;

            ::calculate_exchange_fields(0,::atoms::num_atoms);

            std::cout << cu::x_total_spin_field_array[0] << "\t\t"
                     <<  cu::y_total_spin_field_array[0] << "\t\t"
                     <<  cu::z_total_spin_field_array[0] << "\t\t"
                     << ::atoms::x_total_spin_field_array[0] << "\t\t"
                     << ::atoms::y_total_spin_field_array[0] << "\t\t"
                     << ::atoms::z_total_spin_field_array[0] << "\t\t"
                     << std::endl;

            check_cuda_errors(__FILE__,__LINE__);
            return EXIT_SUCCESS;
         }
      } // end namespace exchange

   } // end namespace internal

} // end namespace vcuda

