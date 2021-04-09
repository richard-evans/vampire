//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
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

//------------------------------------------------------------------------
// Class function to normalise exchange interactions
//------------------------------------------------------------------------
void unitcell::exchange_template_t::normalise_exchange(){

   // Select program to run
   switch(uc::internal::exchange_function){

      case internal::nearest_neighbour:
         return;
         break;

      case internal::shell:{
         return;
         break;
      }

      case internal::exponential:
         normalise_exponential_exchange();
         break;

      default:
         return;

   }

   return;

}

//------------------------------------------------------------------------
// Function to normalise exchange interactions
//------------------------------------------------------------------------
void unitcell::exchange_template_t::normalise_exponential_exchange(){

   // calculate expected sum from all nearest neighbours
   double expected_sum = 0.0;
   for(int a=0; a<ni.size(); a++){
      expected_sum += double(ni[a]);
   }

   // calculate actual sum of all interactions
   double sum = 0.0;
   for(int i=0; i<interaction.size(); i++){
      sum += interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // normalise to get same sum as for nearest neighbours
   const double inv_norm_factor = expected_sum/sum;
   for(int i=0; i<interaction.size(); i++){
      interaction[i].Jij[0][0] *= inv_norm_factor;
      interaction[i].Jij[1][1] *= inv_norm_factor;
      interaction[i].Jij[2][2] *= inv_norm_factor;
   }

   double nsum = 0.0;
   for(int i=0; i<interaction.size(); i++){
      nsum += interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // output sum and normalised sum to screen
   //std::cout << expected_sum << "\t" << sum << "\t" << nsum << std::endl;

   return;

}

} // end if namespace unitcell
