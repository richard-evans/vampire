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
#include <algorithm>
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
   const int num_cells_x_sl  = int(st::internal::total_num_cells) * num_sublattices;

   // variable to compute sum of inverse resistances over all stacks
   double sum_inv_resistance = 0.0;

   // declare memory for storing initial sublattice magnetisations
   std::vector<double> sl_mix(num_sublattices);
   std::vector<double> sl_miy(num_sublattices);
   std::vector<double> sl_miz(num_sublattices);

   //---------------------------------------------------------------------------------------------------------
   // Zero spin torque arrray for parallel version to allow reduction
   //---------------------------------------------------------------------------------------------------------
   #ifdef MPICF
      std::fill(st::internal::cell_sl_spin_torque_fields_x.begin(), st::internal::cell_sl_spin_torque_fields_x.end(), 0.0);
      std::fill(st::internal::cell_sl_spin_torque_fields_y.begin(), st::internal::cell_sl_spin_torque_fields_y.end(), 0.0);
      std::fill(st::internal::cell_sl_spin_torque_fields_z.begin(), st::internal::cell_sl_spin_torque_fields_z.end(), 0.0);
      //std::fill(st::internal::cell_resistance.begin(), st::internal::cell_resistance.end(), 0.0);
      //std::fill(st::internal::stack_resistance.begin(), st::internal::stack_resistance.end(), 0.0); // needed for data output only
      //std::fill(st::internal::stack_current.begin(),    st::internal::stack_current.end(),    0.0);
   #endif

   //---------------------------------------------------------------------------------------------------------
   // loop over all stacks to calculate stack resistance (Parallel loop in MPI)
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

         // save cell resistance to calculate local voltage drop over the cell, and thereby calculate sublattice currents
         st::internal::cell_resistance[cell] = cell_resistance;

      }

      //-----------------------------------------------------
      // accumulate total inverse resistance
      //-----------------------------------------------------
      sum_inv_resistance += 1.0 / total_stack_resistance;

      //------------------------------------------------------------------------------------
      // Compute stack current (as this only depends on the resistance in each stack, R_i)
      //------------------------------------------------------------------------------------
      const double I = st::internal::voltage * program::fractional_electric_field_strength / total_stack_resistance;

      //---------------------------------------------------------
      // Compute cell spin torque fields based on stack currents
      //---------------------------------------------------------
      // // old loop for positive current only
      // for(unsigned int cell = start+1 ; cell < end ; cell++){
      //---------------------------------------------------------
      // loop over all other cells in stack
      // slightly unsafe loop structure, but designed for allowing forward and backward loops starting one after the first cell
      for(unsigned int cell = start + cell_inc ; cell != end + cell_inc ; cell += cell_inc){

         // Calculate vocal voltage drop (based on cell resistance)
         const double V_cell = I * st::internal::cell_resistance[cell];

         // loop over all sublattices, multiply by V_cell due to I_sl = V_cell / R_sl,
         // as torque already pre-divided by R_sl
         const int base = int(cell) * num_sublattices;
         for( int sl = 0; sl < num_sublattices; sl++ ){
            st::internal::cell_sl_spin_torque_fields_x[base + sl] *= V_cell;
            st::internal::cell_sl_spin_torque_fields_y[base + sl] *= V_cell;
            st::internal::cell_sl_spin_torque_fields_z[base + sl] *= V_cell;
         }
      }

      //-----------------------------------------------------
      // save stack resistance and current to arrays
      //-----------------------------------------------------
      st::internal::stack_resistance[stack] = total_stack_resistance;
      st::internal::stack_current[stack]    = I;

   } // end of stack loop

   //------------------------------------------------------------------------------------------
   // Reduce cell spin torque fields and stack currents and resistances on all processors
   // Single reduction per component over the entire flat array
   //------------------------------------------------------------------------------------------
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_spin_torque_fields_x[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_spin_torque_fields_y[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &st::internal::cell_sl_spin_torque_fields_z[0], num_cells_x_sl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sum_inv_resistance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

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

      // hbar / ( 2 e mu_B) = 1.054571817e-34 /(2.0 * 1.602176634e-19 * 9.2740100657e-24) = 35486911.9121
      const double hbar_o_2emuB = 35486911.9121;

      // variable to hold the sum of inverse resistances for all sublattices
      double sum_inv_resistances = 0.0;

      // flat base index for this cell: sl is the fastest-varying index
      const int base = cell * num_sublattices;

      // loop over all sublattices
      for( int sl = 0; sl < num_sublattices; sl++ ){

         const int idx = base + sl;

         const double mix = sl_mix[sl];
         const double miy = sl_miy[sl];
         const double miz = sl_miz[sl];

         // get standard and spin-dependent resistances
         const double Rep = st::internal::cell_sl_resistance[idx];
         const double Rsp = st::internal::cell_sl_spin_resistance[idx];

         // if the sublattice is magnetic then add to sum of inverse resistances
         if(st::internal::sl_magnetic[idx]){

            // calculate next cell reduced magnetization
            const double jsat      = st::internal::cell_sl_isaturation[idx];
            const double mjx       = st::internal::cell_sl_magnetization_x[idx] * jsat;
            const double mjy       = st::internal::cell_sl_magnetization_y[idx] * jsat;
            const double mjz       = st::internal::cell_sl_magnetization_z[idx] * jsat;
            const double alpha     = st::internal::cell_sl_alpha[idx];
            const double mi_dot_mj = ( mix*mjx + miy*mjy + miz*mjz );

            // calculate resistance and sum as 1/R since Rep and Rsp are serial resistances
            const double R_sl = Rep + 0.5*Rsp*(1.0 - mi_dot_mj);
            sum_inv_resistances += 1.0 / R_sl;

            // calculate relative contributions of adiabatic and non-adiabatic spin torque
            const double strj = st::internal::cell_sl_relaxation_torque_rj[idx]; // eta parameter
            const double stpj = st::internal::cell_sl_precession_torque_pj[idx]; // beta parameter (should be multiplied by eta?)

            // calculate field without current based on relative magnetization orientations
            const double hx = (strj-alpha*stpj)*(mjy*miz - mjz*miy) + (stpj+alpha*strj)*mix;
            const double hy = (strj-alpha*stpj)*(mjz*mix - mjx*miz) + (stpj+alpha*strj)*miy;
            const double hz = (strj-alpha*stpj)*(mjx*miy - mjy*mix) + (stpj+alpha*strj)*miz;

            // save field (without current factor, but scaled with R_sl, with I = V/R, so multiply by V_cell to get local current)
            st::internal::cell_sl_spin_torque_fields_x[idx] = hx * jsat * hbar_o_2emuB / R_sl; // multiply by inverse moment
            st::internal::cell_sl_spin_torque_fields_y[idx] = hy * jsat * hbar_o_2emuB / R_sl;
            st::internal::cell_sl_spin_torque_fields_z[idx] = hz * jsat * hbar_o_2emuB / R_sl;

            // update magnetization for next iteration
            sl_mix[sl] = mjx;
            sl_miy[sl] = mjy;
            sl_miz[sl] = mjz;

         }
         else{
            // otherwise just add the standard inverse resistances (skip empty sublattices with R=0)
            if(Rep > 0.0) sum_inv_resistances += 1.0 / Rep;

            // here we would also add the spin resistance here for accumulation in the
            // next step to include the effects of TMR, but for the sublattice
            // calculation this is information is lost and its not clear how this impacts
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

      // flat base index for this cell: sl is the fastest-varying index
      const int base = cell * num_sublattices;

      // loop over all sublattices
      for( int sl = 0; sl < num_sublattices; sl++ ){

         const int idx = base + sl;

         // if the sublattice is magnetic then add to sum of inverse resistances
         if(st::internal::sl_magnetic[idx]){

            // calculate next cell reduced magnetization
            const double isat = st::internal::cell_sl_isaturation[idx];
            sl_mix[sl] = st::internal::cell_sl_magnetization_x[idx] * isat;
            sl_miy[sl] = st::internal::cell_sl_magnetization_y[idx] * isat;
            sl_miz[sl] = st::internal::cell_sl_magnetization_z[idx] * isat;

         }

         // skip empty sublattices (R=0) to avoid short-circuiting the parallel network
         if(st::internal::cell_sl_resistance[idx] > 0.0)
            sum_inv_resistances += 1.0 / st::internal::cell_sl_resistance[idx];

      } // end of sublattice loop

      // calculate total resistance for the cell
      const double R = 1.0 / sum_inv_resistances;

      return R;

   }

} // end of internal namespace
} // end of spin_transport namespace
