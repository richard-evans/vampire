#ifndef __SPIN_FIELDS_HPP__
#define __SPIN_FIELDS_HPP__

// Local cuda headers
#include "cuda_utils.hpp"
#include "exchange_fields.hpp"
#include "data.hpp"
#include "internal.hpp"


namespace vcuda{

    namespace internal{

        __device__ cu_real_t uniaxial_anisotropy_energy(cu::material_parameters_t &material,
                cu_real_t sx, cu_real_t sy, cu_real_t sz);

    }
}


#endif
