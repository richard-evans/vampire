#ifndef __CUDA_EXCHANGE_HPP__
#define __CUDA_EXCHANGE_HPP__

#include "internal.hpp"

#include <iostream>

namespace vcuda
{
   namespace internal
   {
      namespace exchange
      {

          extern int *d_csr_rows;
          extern int *d_coo_rows;
          extern int *d_coo_cols;
          extern cu_real_t *d_coo_vals;


         int initialise_exchange();

         int finalise_exchange();

         int calculate_exchange_fields();


      } // end namespace exchange
   } // end namespace internal
} // end namespace vcuda

#endif
