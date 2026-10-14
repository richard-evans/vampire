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
#include "program.hpp"
#include "spintransport.hpp"
#include "vmpi.hpp"

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
   // Zero spin torque arrray for parallel version to allow reduction
   //---------------------------------------------------------------------------------------------------------
   #ifdef MPICF
      std::fill(spin_transport::internal::cell_spin_torque_fields.begin(), spin_transport::internal::cell_spin_torque_fields.end(), 0.0);
      //std::fill(spin_transport::internal::stack_resistance.begin(), spin_transport::internal::stack_resistance.end(), 0.0); // needed for data output only
      //std::fill(spin_transport::internal::stack_current.begin(),    spin_transport::internal::stack_current.end(),    0.0);
   #endif

   //---------------------------------------------------------------------------------------------------------
   // loop over all stacks to calculate stack resistance (parallelised for MPI)
   //---------------------------------------------------------------------------------------------------------
   for(uint64_t stack = spin_transport::internal::first_stack; stack < spin_transport::internal::last_stack; stack++){

      const unsigned int start = stack_start_index[stack];
      const unsigned int end   = stack_final_index[stack];
      const double isat = spin_transport::internal::cell_isaturation[start]; // saturation magnetization for cell i

      // hbar / ( 2 e mu_B) = 1.054571817e-34 /(2.0 * 1.602176634e-19 * 9.2740100657e-24) = 35486911.9121
      const double hbar_o_2emuB = 35486911.9121;

      double total_stack_resistance = 0.0;

      // load first cell reduced magnetization
      double mix = spin_transport::internal::cell_magnetization[3*start+0] * isat;
      double miy = spin_transport::internal::cell_magnetization[3*start+1] * isat;
      double miz = spin_transport::internal::cell_magnetization[3*start+2] * isat;

      // load first cell resistances for P and AP states
      double Rep = spin_transport::internal::cell_resistance[start];      // electron-phonon scattering resistance
      double Rsp = spin_transport::internal::cell_spin_resistance[start]; // electron-spin scattering resistance

      //------------------------------------------------------------------------------------------------------
      // loop over all other cells in stack starting at cell start+1
      //------------------------------------------------------------------------------------------------------
      //for(unsigned int cell = start+1 ; cell < end ; cell++){
      const int cell_inc = spin_transport::internal::cell_increment;
      // slightly unsafe loop structure, but designed for allowing forward and backward loops starting one after the first cell
      for(unsigned int cell = start + cell_inc ; cell != end + cell_inc ; cell += cell_inc){

         if(spin_transport::internal::magnetic[cell]){
            // calculate next cell reduced magnetization
            const double jsat = spin_transport::internal::cell_isaturation[cell];
            const double mjx = spin_transport::internal::cell_magnetization[3*cell+0] * jsat;
            const double mjy = spin_transport::internal::cell_magnetization[3*cell+1] * jsat;
            const double mjz = spin_transport::internal::cell_magnetization[3*cell+2] * jsat;
            const double alpha = spin_transport::internal::cell_alpha[cell];
            const double mi_dot_mj = ( mix*mjx + miy*mjy + miz*mjz );

            // calculate resistance (need to include T dependence of Rep here)
            total_stack_resistance += Rep + 0.5*Rsp*(1.0 - mi_dot_mj);

            // calculate relavtive contributions of adiabatic and non-adiabatic spin torque
            const double strj = spin_transport::internal::cell_relaxation_torque_rj[cell];
            const double stpj = spin_transport::internal::cell_precession_torque_pj[cell];

            // calculate field without current based on relative magnetization orientations
            const double hx = (strj-alpha*stpj)*(mjy*miz - mjz*miy) + (stpj+alpha*strj)*mix;
            const double hy = (strj-alpha*stpj)*(mjz*mix - mjx*miz) + (stpj+alpha*strj)*miy;
            const double hz = (strj-alpha*stpj)*(mjx*miy - mjy*mix) + (stpj+alpha*strj)*miz;

            // save field (without current factor)
            spin_transport::internal::cell_spin_torque_fields[3*cell+0] = hx * jsat * hbar_o_2emuB; // multiply by inverse moment
            spin_transport::internal::cell_spin_torque_fields[3*cell+1] = hy * jsat * hbar_o_2emuB;
            spin_transport::internal::cell_spin_torque_fields[3*cell+2] = hz * jsat * hbar_o_2emuB;

            // update cell resistances and magnetization
            mix = mjx;
            miy = mjy;
            miz = mjz;

            Rep = spin_transport::internal::cell_resistance[cell];
            Rsp = spin_transport::internal::cell_spin_resistance[cell];

         }
         else{

            // update cell resistances and defer spin resistance to later by accumulating
            Rep += spin_transport::internal::cell_resistance[cell];
            Rsp += spin_transport::internal::cell_spin_resistance[cell]; // The += is significant!

         }

      }

      //-----------------------------------------------------
      // add final non-spin resistance to stack
      //-----------------------------------------------------
      total_stack_resistance += Rep;

      // accumulate total inverse resistance
      sum_inv_resistance += 1.0 / total_stack_resistance;

      //-----------------------------------------------------
      // Compute stack current
      //-----------------------------------------------------
      const double I = spin_transport::internal::voltage * program::fractional_electric_field_strength / total_stack_resistance;

      //---------------------------------------------------------
      // Compute cell spin torque fields based on stack currents
      //---------------------------------------------------------
      // // old loop for positive current only
      // for(unsigned int cell = start+1 ; cell < end ; cell++){
      //---------------------------------------------------------
      // loop over all other cells in stack
      // slightly unsafe loop structure, but designed for allowing forward and backward loops starting one after the first cell
      for(unsigned int cell = start + cell_inc ; cell != end + cell_inc ; cell += cell_inc){
         spin_transport::internal::cell_spin_torque_fields[3*cell+0] *= I;
         spin_transport::internal::cell_spin_torque_fields[3*cell+1] *= I;
         spin_transport::internal::cell_spin_torque_fields[3*cell+2] *= I;
      }

      //-----------------------------------------------------
      // save stack resistance and current to arrays
      //-----------------------------------------------------
      spin_transport::internal::stack_resistance[stack] = total_stack_resistance;
      spin_transport::internal::stack_current[stack]    = I;

   } // end of stack loop

   //------------------------------------------------------------------------------------------
   // Reduce cell spin torque fields and stack currents and resistances on all processors
   //------------------------------------------------------------------------------------------
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &spin_transport::internal::cell_spin_torque_fields[0], 3*spin_transport::internal::total_num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sum_inv_resistance,                       1,                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // save total resistance and current
   spin_transport::total_resistance = 1.0 / sum_inv_resistance;
   spin_transport::total_current = program::fractional_electric_field_strength * spin_transport::internal::voltage / spin_transport::total_resistance;

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
