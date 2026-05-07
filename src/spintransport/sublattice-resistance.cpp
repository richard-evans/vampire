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
#include "program.hpp"
#include "spintransport.hpp"
#include "vmpi.hpp"

// spintransport module headers
#include "internal.hpp"

namespace spin_transport{
namespace internal{

   double calculate_start_cell_resistance_and_magnetization(const int cell,
                                                            std::vector<double>& sl_mix, // sublattice moments of previous cell
                                                            std::vector<double>& sl_miy,
                                                            std::vector<double>& sl_miz);

   // function to calculate cell resistance and torques
   double calculate_cell_resistance_and_torques(const int cell,
                                                std::vector<double>& sl_mix, // sublattice moments of previous cell
                                                std::vector<double>& sl_miy,
                                                std::vector<double>& sl_miz);

//---------------------------------------------------------------------------------------------------------
// Function to calculate stack resistances
//---------------------------------------------------------------------------------------------------------
void calculate_sublattice_resistance(){

   // define local constants for number of cells sublattices to avoid repeated access to global variables
   const int num_sublattices = st::internal::num_sublattices;
   //const int num_cells = st::internal::total_num_cells;

   // variable to compute sum of inverse resistances over all stacks
   double sum_inv_resistance = 0.0;

   // declare memory for storing initial sublattice magnetisations
   std::vector<double> sl_mix(num_sublattices);
   std::vector<double> sl_miy(num_sublattices);
   std::vector<double> sl_miz(num_sublattices);

   //---------------------------------------------------------------------------------------------------------
   // Zero spin torque arrray for parallel version to allow reduction
   //---------------------------------------------------------------------------------------------------------
   //#ifdef MPICF
   //   std::fill(st::internal::cell_spin_torque_fields.begin(), st::internal::cell_spin_torque_fields.end(), 0.0);
      //std::fill(st::internal::stack_resistance.begin(), st::internal::stack_resistance.end(), 0.0); // needed for data output only
      //std::fill(st::internal::stack_current.begin(),    st::internal::stack_current.end(),    0.0);
   //#endif

   // TODO need to parallelise stack loop
   //---------------------------------------------------------------------------------------------------------
   // loop over all stacks to calculate stack resistance (can OpenMP this loop)
   //---------------------------------------------------------------------------------------------------------
   for(uint64_t stack = st::internal::first_stack; stack < st::internal::last_stack; stack++){

      const unsigned int start = stack_start_index[stack];
      const unsigned int end   = stack_final_index[stack];

      // initialise resistance with the first cell
      double total_stack_resistance = calculate_start_cell_resistance_and_magnetization(start, sl_mix, sl_miy, sl_miz);

      //------------------------------------------------------------------------------------------------------
      // loop over all other cells in stack starting at cell start+1
      //------------------------------------------------------------------------------------------------------
      const int cell_inc = st::internal::cell_increment;

      // slightly unsafe loop structure, but designed for allowing forward and backward loops starting one after the first cell
      for(unsigned int cell = start + cell_inc ; cell != end + cell_inc ; cell += cell_inc){

         // calculate resistance of cell and effective torques on each sublattice
         double cell_resistance = calculate_cell_resistance_and_torques(cell, sl_mix, sl_miy, sl_miz);

         // add cell resistance to total resistance
         total_stack_resistance += cell_resistance;

      }

      //-----------------------------------------------------
      // accumulate total inverse resistance
      //-----------------------------------------------------
      sum_inv_resistance += 1.0 / total_stack_resistance;

      //------------------------------------------------------------------------------------
      // Compute stack current (as this only depends on the resistance in each stack, R_i)
      //------------------------------------------------------------------------------------
      const double je = st::internal::voltage * program::fractional_electric_field_strength / total_stack_resistance;

      //---------------------------------------------------------
      // Compute cell spin torque fields based on stack currents
      //---------------------------------------------------------
      // // old loop for positive current only
      // for(unsigned int cell = start+1 ; cell < end ; cell++){
      //---------------------------------------------------------
      // loop over all other cells in stack
      // slightly unsafe loop structure, but designed for allowing forward and backward loops starting one after the first cell
      for(unsigned int cell = start + cell_inc ; cell != end + cell_inc ; cell += cell_inc){
         // loop over all sublattices
         for( int sl = 0; sl < num_sublattices; sl++ ){
            st::internal::cell_sl_spin_torque_fields_x[cell][sl] *= je;
            st::internal::cell_sl_spin_torque_fields_y[cell][sl] *= je;
            st::internal::cell_sl_spin_torque_fields_z[cell][sl] *= je;
         }
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
   // not needed as stack not parallelised and all source data is the same everywhere
   //#ifdef MPICF
   //   MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_spin_torque_fields[0], 3*st::internal::total_num_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   //   //MPI_Allreduce(MPI_IN_PLACE, &st::internal::stack_resistance[0],        st::internal::num_stacks,        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   //   MPI_Allreduce(MPI_IN_PLACE, &sum_inv_resistance,                       1,                               MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   //#endif

   // save total resistance and current
   st::total_resistance = 1.0 / sum_inv_resistance;
   st::total_current = program::fractional_electric_field_strength * st::internal::voltage / st::total_resistance;

   return;

}

   //-----------------------------------------------------------------------------------
   // Function to calculate the resistance of a single cell, forming a set of parallel
   // resistors
   //-----------------------------------------------------------------------------------
   double calculate_cell_resistance_and_torques(const int cell,
                                                std::vector<double>& sl_mix, // sublattice moments of previous cell
                                                std::vector<double>& sl_miy,
                                                std::vector<double>& sl_miz){

      // define local constants for number of cells sublattices to avoid repeated access to global variables
      const int num_sublattices = st::internal::num_sublattices;

      // variable to hold the sum of inverse resistances for all sublattices
      double sum_inv_resistances = 0.0;

      // loop over all sublattices
      for( int sl = 0; sl < num_sublattices; sl++ ){

         const double mix = sl_mix[sl];
         const double miy = sl_miy[sl];
         const double miz = sl_miz[sl];

         // get standard and spin-dependent resistances
         const double Rep = st::internal::cell_sl_resistance[cell][sl];
         const double Rsp = st::internal::cell_sl_spin_resistance[cell][sl];

         // if the sublattice is magnetic then add to sum of inverse resistances
         if(st::internal::sl_magnetic[cell][sl]){

            // calculate next cell reduced magnetization
            const double jsat = st::internal::cell_sl_isaturation[cell][sl];
            const double mjx = st::internal::cell_sl_magnetization_x[cell][sl] * jsat;
            const double mjy = st::internal::cell_sl_magnetization_y[cell][sl] * jsat;
            const double mjz = st::internal::cell_sl_magnetization_z[cell][sl] * jsat;
            const double alpha = st::internal::cell_sl_alpha[cell][sl];
            const double mi_dot_mj = ( mix*mjx + miy*mjy + miz*mjz );

            // calculate resistance and sum as 1/R since Rep and Rsp are serial resistances
            sum_inv_resistances += 1.0 / ( Rep + 0.5*Rsp*(1.0 - mi_dot_mj) );

            // calculate relavtive contributions of adiabatic and non-adiabatic spin torque
            const double strj = st::internal::cell_sl_relaxation_torque_rj[cell][sl];
            const double stpj = st::internal::cell_sl_precession_torque_pj[cell][sl];

            // calculate field without current based on relative magnetization orientations
            const double hx = (strj-alpha*stpj)*(mjy*miz - mjz*miy) + (stpj+alpha*strj)*mix;
            const double hy = (strj-alpha*stpj)*(mjz*mix - mjx*miz) + (stpj+alpha*strj)*miy;
            const double hz = (strj-alpha*stpj)*(mjx*miy - mjy*mix) + (stpj+alpha*strj)*miz;

            // save field (without current factor)
            st::internal::cell_sl_spin_torque_fields_x[cell][sl] = hx;
            st::internal::cell_sl_spin_torque_fields_y[cell][sl] = hy;
            st::internal::cell_sl_spin_torque_fields_z[cell][sl] = hz;

            // update magnetization for next iteration
            sl_mix[sl] = mjx;
            sl_miy[sl] = mjy;
            sl_miz[sl] = mjz;

         }
         else{
            // otherwise just add the standard inverse resistances
            sum_inv_resistances += 1.0 / Rep;

            // here we would also add the spin resistance here for accumulation in the
            // next step to include the effects of TMR, but for the sublattice
            // calcultaion this is information is lost and its not clear how this impacts
            // the cell level resistances that must be done with the inverse sum
            // sum_inv_spin_resistances += 1.0 / Rsp;

         }
      } // end of sublattice loop

      // calculate total resistance for the cell
      const double total_inverse_resistance = sum_inv_resistances;
      const double R = 1.0 / total_inverse_resistance;

      return R;

   }

   //-----------------------------------------------------------------------------------
   // Function to calculate the resistance of the first cell in the stack
   // Each sublattice forms a set of parallel resistors
   //-----------------------------------------------------------------------------------
   double calculate_start_cell_resistance_and_magnetization(const int cell,
                                                            std::vector<double>& sl_mix, // sublattice moments of previous cell
                                                            std::vector<double>& sl_miy,
                                                            std::vector<double>& sl_miz){

      // define local constants for number of cells sublattices to avoid repeated access to global variables
      const int num_sublattices = st::internal::num_sublattices;

      // variable to hold the sum of inverse resistances for all sublattices
      double sum_inv_resistances = 0.0;

      // loop over all sublattices
      for( int sl = 0; sl < num_sublattices; sl++ ){

         // if the sublattice is magnetic then add to sum of inverse resistances
         if(st::internal::sl_magnetic[cell][sl]){

            // calculate next cell reduced magnetization
            const double isat = st::internal::cell_sl_isaturation[cell][sl];
            sl_mix[sl] = st::internal::cell_sl_magnetization_x[cell][sl] * isat;
            sl_miy[sl] = st::internal::cell_sl_magnetization_y[cell][sl] * isat;
            sl_miz[sl] = st::internal::cell_sl_magnetization_z[cell][sl] * isat;

         }

         // otherwise just add the standard inverse resistances
         sum_inv_resistances += 1.0 / st::internal::cell_sl_resistance[cell][sl];

      } // end of sublattice loop

      // calculate total resistance for the cell
      const double R = 1.0 / sum_inv_resistances;

      return R;

   }

} // end of internal namespace
} // end of spin_transport namespace
