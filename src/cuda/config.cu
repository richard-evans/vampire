//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "cuda.hpp"
#include "typedefs.hpp"


// Local cuda headers
#include "cuda_utils.hpp"
#include "data.hpp"

#ifdef CUDA
namespace vcuda{

   namespace config{

      void synchronise(){

         // copy spin data to CPU (via the transfer buffers, which also converts single precision spins)
         vcuda::transfer_spin_positions_from_gpu_to_cpu();

         // configuration files include the magnetisation statistics
         vcuda::stats::get();

         return;

      }

   } // end of namespace config

} // end of namespace vcuda
#endif
