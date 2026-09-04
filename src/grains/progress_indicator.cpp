//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "grains.hpp"
#include "vio.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{

   namespace internal{

      //--------------------------------------------------------------------
      // Prints a "." to stdout and the log at roughly ten equally spaced
      // points while looping over total_grains grains. Guards against the
      // division by zero seen for fewer than ten grains.
      //--------------------------------------------------------------------
      void print_grain_progress(unsigned int grain, unsigned int total_grains){

         const unsigned int stride = total_grains / 10;

         if(stride == 0 || (grain % stride) == 0){
            std::cout << "." << std::flush;
            zlog << "." << std::flush;
         }

         return;

      }

   } // end of internal namespace

} // end of grains namespace
