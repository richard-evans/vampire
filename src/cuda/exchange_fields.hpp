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

      } // end namespace exchange
   } // end namespace internal
} // end namespace vcuda

#endif
