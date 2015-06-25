#include "exchange_fields.hpp"

#include "cuda_utils.hpp"

#include <vector>
#include <thrust/device_vector.h>

namespace atoms
{
   int num_atoms = 100000;
   int total_num_neighbours = 9;

   int exchange_type = 1;

   std::vector<int> nbr_inds;
   std::vector<int> nbr_list;
   std::vector<double> J_vals;

}

namespace vcuda
{

   namespace internal
   {

      bool exchange_initialised = false;

      cusparseHandle_t  cusparse_handle;

      // Diagonal components
      cusparseHybMat_t  Jxx;
      cusparseHybMat_t  Jyy;
      cusparseHybMat_t  Jzz;
      bool Jxx_initialised = false;
      bool Jyy_initialised = false;
      bool Jzz_initialised = false;

      cusparseMatDescr_t Jxx_descr;
      cusparseMatDescr_t Jyy_descr;
      cusparseMatDescr_t Jzz_descr;


      // Off-diagonal components
      cusparseHybMat_t  Jxy;
      cusparseHybMat_t  Jxz;
      cusparseHybMat_t  Jyz;


      int initialise_exchange()
      {

         // Initialise the cuSparse handle
         cusparse_call( cusparseCreate( &cusparse_handle) );

         thrust::device_vector<int> rowptr_d( atoms::nbr_inds);
         thrust::device_vector<int> colinds_d( atoms::nbr_list);
         thrust::device_vector<double> Jvals_d( atoms::J_vals);
         thrust::device_vector<double> Jzvals_d( atoms::J_vals);



         int* rowptr_dptr = thrust::raw_pointer_cast( &rowptr_d[0]);
         int* colinds_dptr = thrust::raw_pointer_cast( &colinds_d[0]);
         double* Jvals_dptr = thrust::raw_pointer_cast( &Jvals_d[0]);
         double* Jzvals_dptr = thrust::raw_pointer_cast( &Jzvals_d[0]);
         cusparseStatus_t status;

         size_t avail;
         size_t total;
         cudaMemGetInfo( &avail, &total );
         size_t used = total - avail;
         std::cout << "Total: " << total << " , Used: " << used << " , Free: " << avail << std::endl;

         switch( atoms::exchange_type)
         {
            case 0: // Isotropic

               //--------------------------------------------------------------
               // Exchange is isotropic so Jxx = Jyy = Jzz
               // and Jxy = Jxz = Jyx = 0
               //--------------------------------------------------------------
               check_cuda_errors(__FILE__,__LINE__);


               // Initialise the  hybrid matrix
               cusparseCreateHybMat( &Jxx);


               cusparseCreateMatDescr(&Jxx_descr);
               cusparseSetMatType(Jxx_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
               cusparseSetMatIndexBase(Jxx_descr,CUSPARSE_INDEX_BASE_ZERO);

               status = cusparseDcsr2hyb(
                     cusparse_handle,
                     atoms::num_atoms,
                     atoms::num_atoms,
                     Jxx_descr,
                     Jvals_dptr,
                     rowptr_dptr,
                     colinds_dptr,
                     Jxx,
                     atoms::total_num_neighbours/ atoms::num_atoms,
                     CUSPARSE_HYB_PARTITION_AUTO);


               if( status != CUSPARSE_STATUS_SUCCESS)
               {
                  std::cerr << "Error: Failed to conver csr to hyb format" << std::endl;
               }

               Jxx_initialised = true;
               check_cuda_errors(__FILE__,__LINE__);

               break;

            case 1: // Vector


               cusparseCreateHybMat( &Jxx);

               cusparseCreateMatDescr(&Jxx_descr);
               cusparseSetMatType(Jxx_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
               cusparseSetMatIndexBase(Jxx_descr,CUSPARSE_INDEX_BASE_ZERO);

               status = cusparseDcsr2hyb(
                     cusparse_handle,
                     atoms::num_atoms,
                     atoms::num_atoms,
                     Jxx_descr,
                     Jvals_dptr,
                     rowptr_dptr,
                     colinds_dptr,
                     Jxx,
                     atoms::total_num_neighbours/ atoms::num_atoms,
                     CUSPARSE_HYB_PARTITION_AUTO);


               if( status != CUSPARSE_STATUS_SUCCESS)
               {
                  std::cerr << "Error: Failed to conver csr to hyb format" << std::endl;
               }

               Jxx_initialised = true;


               cusparseCreateHybMat( &Jyy);

               cusparseCreateMatDescr(&Jyy_descr);
               cusparseSetMatType(Jyy_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
               cusparseSetMatIndexBase(Jyy_descr,CUSPARSE_INDEX_BASE_ZERO);

               status = cusparseDcsr2hyb(
                     cusparse_handle,
                     atoms::num_atoms,
                     atoms::num_atoms,
                     Jyy_descr,
                     Jvals_dptr,
                     rowptr_dptr,
                     colinds_dptr,
                     Jyy,
                     atoms::total_num_neighbours/ atoms::num_atoms,
                     CUSPARSE_HYB_PARTITION_AUTO);


               if( status != CUSPARSE_STATUS_SUCCESS)
               {
                  std::cerr << "Error: Failed to conver csr to hyb format" << std::endl;
               }

               Jyy_initialised = true;

               cusparseCreateHybMat( &Jzz);

               cusparseCreateMatDescr(&Jzz_descr);
               cusparseSetMatType(Jzz_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
               cusparseSetMatIndexBase(Jzz_descr,CUSPARSE_INDEX_BASE_ZERO);

               status = cusparseDcsr2hyb(
                     cusparse_handle,
                     atoms::num_atoms,
                     atoms::num_atoms,
                     Jzz_descr,
                     Jzvals_dptr,
                     rowptr_dptr,
                     colinds_dptr,
                     Jzz,
                     atoms::total_num_neighbours/ atoms::num_atoms,
                     CUSPARSE_HYB_PARTITION_AUTO);


               if( status != CUSPARSE_STATUS_SUCCESS)
               {
                  std::cerr << "Error: Failed to conver csr to hyb format" << std::endl;
               }

               Jzz_initialised = true;

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

            if( Jxx_initialised) cusparseDestroyHybMat( Jxx );
            if( Jyy_initialised) cusparseDestroyHybMat( Jyy );
            if( Jzz_initialised) cusparseDestroyHybMat( Jzz );

         }
         return EXIT_SUCCESS;
      }

      int calculate_exchange_fields()
      {

         // check calling of routine if error checking is activated
         if(true){std::cout << "calculate_exchange_fields has been called. N_atoms = "<< atoms::num_atoms  << std::endl;}

         check_cuda_errors(__FILE__,__LINE__);

         thrust::device_vector<double> Sx;
         thrust::device_vector<double> Sy;
         thrust::device_vector<double> Sz;

         thrust::device_vector<double> Hx;
         thrust::device_vector<double> Hy;
         thrust::device_vector<double> Hz;

         try{

            Sx.assign( atoms::num_atoms, 1.0);
            Sy.assign( atoms::num_atoms, 1.0);
            Sz.assign( atoms::num_atoms, 1.0);

            Hx.assign( atoms::num_atoms, 0.0);
            Hy.assign( atoms::num_atoms, 0.0);
            Hz.assign( atoms::num_atoms, 0.0);
         }
         catch( std::bad_alloc)
         {
            std::cerr << "Bad alloc thrown." << std::endl;
            exit(-1);
         }
         catch( thrust::system_error err)
         {
            std::cerr << "Thrust system error caught: " << err.what() << std::endl;
            exit(-1);
         }

         double alpha_d = 1.0;
         double beta_d = 1.0;

         double* Sx_dptr = thrust::raw_pointer_cast( &Sx[0]);
         double* Sy_dptr = thrust::raw_pointer_cast( &Sy[0]);
         double* Sz_dptr = thrust::raw_pointer_cast( &Sz[0]);

         double* Hx_dptr = thrust::raw_pointer_cast( &Hx[0]);
         double* Hy_dptr = thrust::raw_pointer_cast( &Hy[0]);
         double* Hz_dptr = thrust::raw_pointer_cast( &Hz[0]);

         check_cuda_errors(__FILE__,__LINE__);
         switch( atoms::exchange_type)
         {
            case 0: // Isotropic

               std::cout << "Isotropic Exchange" << std::endl;

               //--------------------------------------------------------------
               // Exchange is isotropic so Jxx = Jyy = Jzz
               // and Jxy = Jxz = Jyx = 0
               //--------------------------------------------------------------

               if( !exchange_initialised) initialise_exchange();

               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jxx_descr,
                     Jxx,
                     Sx_dptr,
                     &beta_d,
                     Hx_dptr);

               check_cuda_errors(__FILE__,__LINE__);
               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jxx_descr,
                     Jxx,
                     Sy_dptr,
                     &beta_d,
                     Hy_dptr);

         check_cuda_errors(__FILE__,__LINE__);
               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jxx_descr,
                     Jxx,
                     Sz_dptr,
                     &beta_d,
                     Hz_dptr);
         check_cuda_errors(__FILE__,__LINE__);
               break;

            case 1: // Vector exchange

               std::cout << "Vector Exchange" << std::endl;
               //--------------------------------------------------------------
               // Exchange is diagonal so Jxx != Jyy != Jzz
               // and Jxy = Jxz = Jyx = 0
               //--------------------------------------------------------------

               if( !exchange_initialised) initialise_exchange();

               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jxx_descr,
                     Jxx,
                     Sx_dptr,
                     &beta_d,
                     Hx_dptr);

               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jyy_descr,
                     Jyy,
                     Sy_dptr,
                     &beta_d,
                     Hy_dptr);

               cusparseDhybmv(
                     cusparse_handle,
                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha_d,
                     Jzz_descr,
                     Jzz,
                     Sz_dptr,
                     &beta_d,
                     Hz_dptr);


               break;

         }
         check_cuda_errors(__FILE__,__LINE__);

         for( int i = 0; i < atoms::num_atoms; i++) std::cerr << i << "\t" << Hx[i] << "\t" << Hy[i] << "\t" << Hz[i] << std::endl;

         std::cout << "Finished field calculation." << std::endl;

         check_cuda_errors(__FILE__,__LINE__);
         return EXIT_SUCCESS;
      }


   } // end namespace internal

} // end namespace vcuda

