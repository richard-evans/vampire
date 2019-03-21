//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "spintransport.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{
namespace internal{

//---------------------------------------------------------------------------------------------------------
// Function to calculate stack resistances
//---------------------------------------------------------------------------------------------------------
void calculate_magnetoresistance(){

   // variable to compute sum of inverse resistances
   double sum_inv_resistance = 0.0;

   //---------------------------------------------------------------------------------------------------------
   // loop over all stacks to calculate stack resistance
   //---------------------------------------------------------------------------------------------------------
   for(int stack = 0; stack < st::internal::num_stacks; stack++){

      const unsigned int start = stack_start_index[stack];
      const unsigned int end   = stack_final_index[stack];

      double total_stack_resistance = 0.0;

      // load first cell magnetization
      double mix = st::internal::cell_magnetization[3*start+0];
      double miy = st::internal::cell_magnetization[3*start+1];
      double miz = st::internal::cell_magnetization[3*start+2];

      // load first cell resistances for P and AP states
      double Rep = st::internal::cell_resistance[start];      // electron-phonon scattering resistance
      double Rsp = st::internal::cell_spin_resistance[start]; // electron-spin scattering resistance

      //------------------------------------------------------------------------------------------------------
      // loop over all other cells in stack starting at cell start+1
      //------------------------------------------------------------------------------------------------------
      for(unsigned int cell = start+1 ; cell < end ; cell++){

         if(st::internal::magnetic[cell]){
            // calculate next cell magnetization
            const double mjx = st::internal::cell_magnetization[3*cell+0];
            const double mjy = st::internal::cell_magnetization[3*cell+1];
            const double mjz = st::internal::cell_magnetization[3*cell+2];

            const double mi_dot_mj = ( mix*mjx + miy*mjy + miz*mjz );

            // calculate resistance (need to include T dependence of Rep here)
            total_stack_resistance += Rep + 0.5*Rsp*(1.0 - mi_dot_mj);

            // update cell resistances and magnetization
            mix = mjx;
            miy = mjy;
            miz = mjz;

            Rep = st::internal::cell_resistance[cell];
            Rsp = st::internal::cell_spin_resistance[cell];

         }
         else{

            // update cell resistances and defer spin resistance to later by accumulating
            Rep += st::internal::cell_resistance[cell];
            Rsp += st::internal::cell_spin_resistance[cell]; // The += is significant!

         }

      }

      //-----------------------------------------------------
      // add final non-spin resistance to stack
      //-----------------------------------------------------
      total_stack_resistance += Rep;

      //-----------------------------------------------------
      // save stack resistance to array
      //-----------------------------------------------------
      st::internal::stack_resistance[stack] = total_stack_resistance;

      sum_inv_resistance += 1.0 / total_stack_resistance;

   } // end of stack loop

   //-----------------------------------------------------
   // Compute stack currents
   //-----------------------------------------------------
   for(int stack = 0; stack < st::internal::num_stacks; stack++){
      st::internal::stack_current[stack] = st::internal::voltage / st::internal::stack_resistance[stack];
   }

   // save total resistance and current
   st::total_resistance = 1.0 / sum_inv_resistance;
   st::total_current = st::internal::voltage / st::total_resistance;

   std::cout << "R: " << st::total_resistance << "\t" << "I: " << st::total_current << std::endl;

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
