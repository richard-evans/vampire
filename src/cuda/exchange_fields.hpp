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



         int initialise_exchange();

         int finalise_exchange();

         int calculate_exchange_fields();


      } // end namespace exchange
   } // end namespace internal
} // end namespace vcuda

#endif
