#ifndef __THERMAL_FIELDS_CUDA_HPP__
#define __THERMAL_FIELDS_CUDA_HPP__

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
