//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016, Jack B Collings 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

//------------------------------------------------------------------------
// Function to determine exchange energy based on interaction range and
// exchange interaction function
//------------------------------------------------------------------------

double exchange(double range_sq, double nn_cutoff_sq){

   // Select program to run
   switch(uc::internal::exchange_function){

      case nearest_neighbour:{
         if(range_sq <= nn_cutoff_sq) return 1.0;
         else return 0.0;
         break;
      }

      case shell:{
         return 1.0;
      }

      case exponential:{
         return uc::internal::exchange_multiplier*exp(-sqrt(range_sq)/uc::internal::exchange_decay) + uc::internal::exchange_shift;
         break;
      }

      default:{
         if(range_sq <= nn_cutoff_sq) return 1.0;
         else return 0.0;
      }
   }

   return 0.0;

}

} // end of internal namespace
} // end of unitcell namespace
