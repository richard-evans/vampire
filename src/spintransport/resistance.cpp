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
      std::fill(st::internal::cell_spin_torque_fields.begin(), st::internal::cell_spin_torque_fields.end(), 0.0);
      //std::fill(st::internal::stack_resistance.begin(), st::internal::stack_resistance.end(), 0.0); // needed for data output only
      //std::fill(st::internal::stack_current.begin(),    st::internal::stack_current.end(),    0.0);
   #endif

   // TODO need to parallelise stack loop
   //---------------------------------------------------------------------------------------------------------
   // loop over all stacks to calculate stack resistance (can OpenMP this loop)
   //---------------------------------------------------------------------------------------------------------
   for(uint64_t stack = st::internal::first_stack; stack < st::internal::last_stack; stack++){

      const unsigned int start = stack_start_index[stack];
      const unsigned int end   = stack_final_index[stack];
      const double isat = st::internal::cell_isaturation[start]; // saturation magnetization for cell i

      double total_stack_resistance = 0.0;

      // load first cell reduced magnetization
      double mix = st::internal::cell_magnetization[3*start+0] * isat;
      double miy = st::internal::cell_magnetization[3*start+1] * isat;
      double miz = st::internal::cell_magnetization[3*start+2] * isat;

      // load first cell resistances for P and AP states
      double Rep = st::internal::cell_resistance[start];      // electron-phonon scattering resistance
      double Rsp = st::internal::cell_spin_resistance[start]; // electron-spin scattering resistance

      //------------------------------------------------------------------------------------------------------
      // loop over all other cells in stack starting at cell start+1
      //------------------------------------------------------------------------------------------------------
      for(unsigned int cell = start+1 ; cell < end ; cell++){

         if(st::internal::magnetic[cell]){
            // calculate next cell reduced magnetization
            const double jsat = st::internal::cell_isaturation[cell];
            const double mjx = st::internal::cell_magnetization[3*cell+0] * jsat;
            const double mjy = st::internal::cell_magnetization[3*cell+1] * jsat;
            const double mjz = st::internal::cell_magnetization[3*cell+2] * jsat;
            const double alpha = st::internal::cell_alpha[cell];
            const double mi_dot_mj = ( mix*mjx + miy*mjy + miz*mjz );

            // calculate resistance (need to include T dependence of Rep here)
            total_stack_resistance += Rep + 0.5*Rsp*(1.0 - mi_dot_mj);

            // calculate relavtive contributions of adiabatic and non-adiabatic spin torque
            const double strj = st::internal::cell_relaxation_torque_rj[cell];
            const double stpj = st::internal::cell_precession_torque_pj[cell];

            // calculate field without current based on relative magnetization orientations
            const double hx = (strj-alpha*stpj)*(mjy*miz - mjz*miy) + (stpj+alpha*strj)*mix;
            const double hy = (strj-alpha*stpj)*(mjz*mix - mjx*miz) + (stpj+alpha*strj)*miy;
            const double hz = (strj-alpha*stpj)*(mjx*miy - mjy*mix) + (stpj+alpha*strj)*miz;

            // save field (without current factor)
            st::internal::cell_spin_torque_fields[3*cell+0] = hx;
            st::internal::cell_spin_torque_fields[3*cell+1] = hy;
            st::internal::cell_spin_torque_fields[3*cell+2] = hz;

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

      // accumulate total inverse resistance
      sum_inv_resistance += 1.0 / total_stack_resistance;

      //-----------------------------------------------------
      // Compute stack current
      //-----------------------------------------------------
      const double je = st::internal::voltage * program::fractional_electric_field_strength / total_stack_resistance;

      //---------------------------------------------------------
      // Compute cell spin torque fields based on stack currents
      //---------------------------------------------------------
      // loop over all other cells in stack starting at cell start+1
      for(unsigned int cell = start+1 ; cell < end ; cell++){
         st::internal::cell_spin_torque_fields[3*cell+0] *= je;
         st::internal::cell_spin_torque_fields[3*cell+1] *= je;
         st::internal::cell_spin_torque_fields[3*cell+2] *= je;
      }

      //-----------------------------------------------------
      // save stack resistance and current to arrays
      //-----------------------------------------------------
      st::internal::stack_resistance[stack] = total_stack_resistance;
      st::internal::stack_current[stack]    = je;

   } // end of stack loop

   //------------------------------------------------------------------------------------------
   // Reduce cell spin trorque fields and stack currents and resistances on all processors
   //------------------------------------------------------------------------------------------
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_spin_torque_fields[0], 3*st::internal::total_num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //MPI_Allreduce(MPI_IN_PLACE, &st::internal::stack_resistance[0],        st::internal::num_stacks,        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sum_inv_resistance,                       1,                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // save total resistance and current
   st::total_resistance = 1.0 / sum_inv_resistance;
   st::total_current = program::fractional_electric_field_strength * st::internal::voltage / st::total_resistance;

   return;

}

} // end of internal namespace
} // end of spin_transport namespace
