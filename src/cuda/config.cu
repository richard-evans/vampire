//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"

// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"

#ifdef CUDA
namespace vcuda{

   namespace config{

      void synchronise(){

		   // copy spin data to CPU
         thrust::copy(internal::atoms::x_spin_array.begin(),internal::atoms::x_spin_array.end(),::atoms::x_spin_array.begin());
         thrust::copy(internal::atoms::y_spin_array.begin(),internal::atoms::y_spin_array.end(),::atoms::y_spin_array.begin());
         thrust::copy(internal::atoms::z_spin_array.begin(),internal::atoms::z_spin_array.end(),::atoms::z_spin_array.begin());

         return;

      }

   } // end of namespace config

} // end of namespace vcuda
#endif
