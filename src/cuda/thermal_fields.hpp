#ifndef CUDA_THERMAL_FIELDS_H_
#define CUDA_THERMAL_FIELDS_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) Andrea Meo 2022. All rights reserved.
//
//------------------------------------------------------------------------------

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"
#include "internal.hpp"


namespace vcuda{

    namespace internal{

      __device__ cu_real_t calculate_thermal_field(const cu_real_t temperature, const cu_real_t alpha,
                                                   const cu_real_t Tc, const cu_real_t sigma,
                                                   curandState local_state);

    }
}


#endif  // CUDA_THERMAL_FIELDS_H_
