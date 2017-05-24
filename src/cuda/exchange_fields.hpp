#ifndef __CUDA_EXCHANGE_HPP__
#define __CUDA_EXCHANGE_HPP__

#include "cusparse.h"

#include <iostream>

namespace vcuda
{
   namespace internal
   {
      namespace exchange
      {

         inline
            void cusparse_call( cusparseStatus_t status)
            {
               if( status != CUSPARSE_STATUS_SUCCESS)
               {
                  std::cerr << "Error: cusparse failed at " << __FILE__ << ", " << __LINE__ << std::endl;
               }
            }


         int initialise_exchange();

         int finalise_exchange();

         int calculate_exchange_fields();

         inline cusparseStatus_t cusparseTcsrmv (
                 cusparseHandle_t handle,
                 cusparseOperation_t transA,
                 int m,
                 int n,
                 int nnz,
                 const double *alpha,
                 const cusparseMatDescr_t descrA,
                 const double *csrValA,
                 const int *csrRowPtrA,
                 const int *csrColIndA,
                 const double *x,
                 const double *beta,
                 double *y
                 )
         {
             return cusparseDcsrmv (
                 handle,
                 transA,
                 m,
                 n,
                 nnz,
                 alpha,
                 descrA,
                 csrValA,
                 csrRowPtrA,
                 csrColIndA,
                 x,
                 beta,
                 y
                 );
         }

         inline cusparseStatus_t cusparseTcsrmv (
                 cusparseHandle_t handle,
                 cusparseOperation_t transA,
                 int m,
                 int n,
                 int nnz,
                 const float *alpha,
                 const cusparseMatDescr_t descrA,
                 const float *csrValA,
                 const int *csrRowPtrA,
                 const int *csrColIndA,
                 const float *x,
                 const float *beta,
                 float *y
                 )
         {
             return cusparseScsrmv (
                 handle,
                 transA,
                 m,
                 n,
                 nnz,
                 alpha,
                 descrA,
                 csrValA,
                 csrRowPtrA,
                 csrColIndA,
                 x,
                 beta,
                 y
                 );
         }

      } // end namespace exchange
   } // end namespace internal
} // end namespace vcuda

#endif
