#ifndef CUDA_SPIN_FIELDS_H_
#define CUDA_SPIN_FIELDS_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) M A Ellis 2021. All rights reserved.
//
//------------------------------------------------------------------------------

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
