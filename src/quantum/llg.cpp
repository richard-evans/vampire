//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Fried-Conrad Weber 2025. All rights reserved.
//
//   Email: fried-conrad.weber@uni-potsdam.de
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "quantum.hpp"
#include "material.hpp"
#include "vio.hpp"

#include "internal.hpp"

#ifdef CUDA
#include "llg-HO-cuda.hpp"
#endif

namespace quantum{

   namespace internal{

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // LLG Wrapper Function
   //---------------------------------------------------------------------------
   void llg(){

      // Call LLG method based on integration type
      #ifdef MPICF
         // MPI parallel version
         if(internal::llg_method == internal::llg_fft){
            internal::llg_FFT_mpi();
         }
         else if(internal::llg_method == internal::llg_ho){
            internal::llg_HO_mpi();
         }
      #else
         // Serial version
         #ifdef CUDA
            // CUDA accelerated version
            if(internal::llg_method == internal::llg_ho){
               internal::cuda_ho::llg_HO_step();
            }
            else if(internal::llg_method == internal::llg_fft){
               internal::llg_FFT();
            }
         #else
            // CPU serial version
            if(internal::llg_method == internal::llg_fft){
               internal::llg_FFT();
            }
            else if(internal::llg_method == internal::llg_ho){
               internal::llg_HO();
            }
         #endif
      #endif

      return;
   }

} // end of quantum namespace
