//------------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) O Arbelaez Echeverri, M A Ellis & R F L Evans 2015. All rights reserved.
//
//------------------------------------------------------------------------------
#ifndef CUDA_MONTE_CARLO_HPP_
#define CUDA_MONTE_CARLO_HPP_

#include "cuda.hpp"
#include "data.hpp"
#include "internal.hpp"
#include "typedefs.hpp"


namespace vcuda {
    namespace internal {
        namespace mc {

            extern bool initialised;

            int initialise();
            void finalise();
            void __mc_step();
        }
    }

} /*vcuda*/

#endif
