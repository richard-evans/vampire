#ifndef CUDA_EXCHANGE_FIELDS_H_
#define CUDA_EXCHANGE_FIELDS_H_
//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   GNU GPL (version 2) licence (see licence file for details).
//
//   (c) M A Ellis 2015. All rights reserved.
//
//------------------------------------------------------------------------------

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

            extern cu_real_t *d_spin3n;

            extern bool exchange_initialised;

            // neighbour list of each atom for isotropic exchange (used instead of the 3N x 3N matrix for spin fields)
            extern bool use_isotropic_kernel;
            extern int *d_nbr_start;
            extern int *d_nbr_atom;
            extern cu_real_t *d_nbr_J;


            int initialise_exchange();

            int finalise_exchange();

            int calculate_exchange_fields();

            __device__ cu_real_t exchange_field_component(int *csr_cols, int* csr_rows, cu_real_t *vals, cu_real_t *spin3n, const int i);

        } // end namespace exchange
    } // end namespace internal
} // end namespace vcuda

#endif
