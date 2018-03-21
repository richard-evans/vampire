//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
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

// forward function declaration
void normalise_exponential_exchange(unitcell::unit_cell_t& unit_cell);

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

      case exponential:{
         return exp(-sqrt(range_sq)/uc::internal::exchange_decay);
         break;
      }

      default:{
         if(range_sq <= nn_cutoff_sq) return 1.0;
         else return 0.0;
      }
   }

   return 0.0;

}

//------------------------------------------------------------------------
// Function to normalise exchange interactions
//------------------------------------------------------------------------
void normalise_exchange(unitcell::unit_cell_t& unit_cell){

   // Select program to run
   switch(uc::internal::exchange_function){

      case nearest_neighbour:
         return;
         break;

      case exponential:
         normalise_exponential_exchange(unit_cell);
         break;

      default:
         return;

   }

   return;

}

//------------------------------------------------------------------------
// Function to normalise exchange interactions
//------------------------------------------------------------------------
void normalise_exponential_exchange(unitcell::unit_cell_t& unit_cell){

   // calculate expected sum from all nearest neighbours
   double expected_sum = 0.0;
   for(int a=0; a<unit_cell.atom.size(); a++){
      expected_sum += double(unit_cell.atom[a].ni);
   }

   // calculate actual sum of all interactions
   double sum = 0.0;
   for(int i=0; i<unit_cell.interaction.size(); i++){
      sum += unit_cell.interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // normalise to get same sum as for nearest neighbours
   const double inv_norm_factor = expected_sum/sum;
   for(int i=0; i<unit_cell.interaction.size(); i++){
      unit_cell.interaction[i].Jij[0][0] *= inv_norm_factor;
      unit_cell.interaction[i].Jij[1][1] *= inv_norm_factor;
      unit_cell.interaction[i].Jij[2][2] *= inv_norm_factor;

      // normalise biquadratic exchange interactions
      unit_cell.biquadratic_interaction[i].Jij[0][0] *= inv_norm_factor;
      unit_cell.biquadratic_interaction[i].Jij[1][1] *= inv_norm_factor;
      unit_cell.biquadratic_interaction[i].Jij[2][2] *= inv_norm_factor;

   }

   double nsum = 0.0;
   for(int i=0; i<unit_cell.interaction.size(); i++){
      nsum += unit_cell.interaction[i].Jij[0][0]; // xx only since J is a trace anyway
   }

   // output sum and normalised sum to screen
   //std::cout << expected_sum << "\t" << sum << "\t" << nsum << std::endl;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
