#ifndef __THERMAL_FIELDS_CUDA_HPP__
#define __THERMAL_FIELDS_CUDA_HPP__
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


#endif  // __THERMAL_FIELDS_CUDA_HPP__
