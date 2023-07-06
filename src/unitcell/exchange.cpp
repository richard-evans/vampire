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

double exchange(double range, double cutoff, int mat_i, int mat_j){

   // Select program to run
   switch(uc::internal::exchange_function){

      case nearest_neighbour:{
         if(range <= cutoff) return 1.0;
         else return 0.0;
         break;
      }

      case shell:{
         return 1.0;
      }

      case exponential:{
         // Set exchange parameters ready for general normalisation function.
         for (size_t i = 0; i < uc::internal::material_exchange_parameters.size(); ++i){
            for (size_t j = 0; j < uc::internal::material_exchange_parameters.size(); ++j){
               if (j >= i){
                  uc::internal::material_exchange_parameters[i][j].decay_multiplier = uc::internal::exchange_multiplier;
                  uc::internal::material_exchange_parameters[i][j].decay_length = uc::internal::exchange_decay;
                  uc::internal::material_exchange_parameters[i][j].decay_shift = uc::internal::exchange_shift;
               }
            }
         }
         return uc::internal::exchange_multiplier*exp(-range/uc::internal::exchange_decay) + uc::internal::exchange_shift;
         break;
      }

      case material_exponential:{
         unsigned int min = std::min(mat_i, mat_j);
         unsigned int max = std::max(mat_i, mat_j);
         double A = uc::internal::material_exchange_parameters[min][max].decay_multiplier; // only need min and max due to symmetry of exchange
         double B = uc::internal::material_exchange_parameters[min][max].decay_length;
         double C = uc::internal::material_exchange_parameters[min][max].decay_shift;
         double Jij = A*exp(-range/B) + C;
         return Jij;
         break;
      }

      case RKKY:{
         return (sin(2*RKKYkf*range) - 2*RKKYkf*range*cos(2*RKKYkf*range))/((RKKYkf*range)*(RKKYkf*range)*(RKKYkf*range)*(RKKYkf*range));
         break;
      }

      default:{
         if(range <= cutoff) return 1.0;
         else return 0.0;
      }
   }

   return 0.0;

}

} // end of internal namespace
} // end of unitcell namespace
