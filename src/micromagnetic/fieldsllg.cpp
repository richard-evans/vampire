//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2021. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "cells.hpp"
#include "../cells/internal.hpp"
#include "constants.hpp"
#include "dipole.hpp"
#include "environment.hpp"
#include "errors.hpp"
#include "micromagnetic.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// micromagnetic module headers
#include "internal.hpp"

// shorthand for brevity
namespace mm = micromagnetic::internal;

//---------------------------------------------------------------------------
// Function to calculate spin-dependent fields for LLG equation
//---------------------------------------------------------------------------
void mm::calculate_llg_spin_fields(const double temperature,
                               const int num_cells,
                               std::vector<double>& mx_array, // Normalised magnetic moment unit vector
                               std::vector<double>& my_array, // Normalised magnetic moment unit vector
                               std::vector<double>& mz_array, // Normalised magnetic moment unit vector
                               std::vector<double>& x_total_spin_field_array,   // total magnetic field
                               std::vector<double>& y_total_spin_field_array,   // total magnetic field
                               std::vector<double>& z_total_spin_field_array){  // total magnetic field

   // Loop over all micromagnetic cells - should this be all cells?
   for (int cell = 0; cell < num_cells; cell++){
      // Optionally determine temperature dependent constants
      if(mm::temperature_dependent_parameters){

         const double rT = temperature/Tc[cell]; // ratio T/Tc

         // determine equilibrium magnetization and damping
         if( temperature <= Tc[cell] ){
            m_e[cell] = std::pow( (1.0 - rT) ,0.34);
            alpha_perp[cell] = alpha[cell]*(1.0 - 0.333333333333 * rT);
         }
         else{
            m_e[cell] = 0.01;
            alpha_perp[cell] = alpha[cell] * 0.666666666667 * rT;
         }

      }
      // constant micromagnetic parameters
      else{
         m_e[cell] = 1.0;
         alpha_perp[cell] = alpha[cell];
      }

   }

   // Determine fields for all micromagnetic cells
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

      // determine cell ID of cell
      int cell = list_of_micromagnetic_cells[lc];

      // temporary constants for brevity
      const double mx = mx_array[cell];
      const double my = my_array[cell];
      const double mz = mz_array[cell];
      const double me = m_e[cell];

      //------------------------------------------------------------------------
      // calculates the exchange fields as me^1.71 *A*(xi-xj)/m_e^2
      //------------------------------------------------------------------------

      // calculate cross-sectional area of cells
      //const double area = cells::macro_cell_size_x * cells::macro_cell_size_y;

      // Temporary variable to calculate total exchange field
      double exchange_field[3]={ 0.0, 0.0, 0.0 };

      // Check that there is more than one cell
      if (num_cells > 1){

         const int start = macro_neighbour_list_start_index[cell]; // save start index for neighbour list
         const int end = macro_neighbour_list_end_index[cell] +1;  // save end index for neighbour list

         // loop over neighbouring cells
         for(int j = start;j< end;j++){

            // get ID of neighbouring cell
            const int cellj = macro_neighbour_list_array[j];

            // get equilibrim magnetization of neighbouring cell
            const double mj = m_e[cellj];

            // calculate reduced exchange constant factor
            double Ac = A[j]*std::pow(mj,1.71);

            // Add field from cell to total exchange field (at equillibrium this term goes to zero)
            exchange_field[0] -= Ac*(mx_array[cellj]*m_e[cellj] - mx*me);
            exchange_field[1] -= Ac*(my_array[cellj]*m_e[cellj] - my*me);
            exchange_field[2] -= Ac*(mz_array[cellj]*m_e[cellj] - mz*me);

         } // end of loop over cells

      } // end of check for neighbouring cells

      // get cell level spin transfer torque parameters
  		const double strj  = mm::stt_rj[cell];
  		const double stpj  = mm::stt_pj[cell];
      const double alpha = alpha_perp[cell]; // get local cell alpha

      // calculate field (temperature?)
		const double hsttx = (strj-alpha*stpj)*(my*mm::sttpz - mz*mm::sttpy) + (stpj+alpha*strj)*mm::sttpx;
		const double hstty = (strj-alpha*stpj)*(mz*mm::sttpx - mx*mm::sttpz) + (stpj+alpha*strj)*mm::sttpy;
		const double hsttz = (strj-alpha*stpj)*(mx*mm::sttpy - my*mm::sttpx) + (stpj+alpha*strj)*mm::sttpz;

      x_total_spin_field_array[cell] = exchange_field[0] - one_o_chi_perp[cell]*mx*me*ku_x[cell] + hsttx;
      y_total_spin_field_array[cell] = exchange_field[1] - one_o_chi_perp[cell]*my*me*ku_y[cell] + hstty;
      z_total_spin_field_array[cell] = exchange_field[2] - one_o_chi_perp[cell]*mz*me*ku_z[cell] + hsttz;

   } // end of loop over cells

   return;

}

void mm::calculate_llg_external_fields(const double temperature,
                                       const int num_cells,
                                       std::vector<double>& x_array,
                                       std::vector<double>& y_array,
                                       std::vector<double>& z_array,
                                       std::vector<double>& x_total_external_field_array,
                                       std::vector<double>& y_total_external_field_array,
                                       std::vector<double>& z_total_external_field_array){

   // temporary constants for brevity
   const double kB = constants::kB;

   // Determine fields for all micromagnetic cells
   for (int lc = 0; lc < number_of_micromagnetic_cells; lc++){

      // determine cell ID of cell
      const int cell = list_of_micromagnetic_cells[lc];

      // calculate width of thermal field
      const double sigma_perp = sqrt( 2.0 * kB * temperature * mm::alpha_perp[cell] / ( mm::ms[cell] * mp::dt ) );

      x_total_external_field_array[cell] = mm::ext_field[0] + sigma_perp*mtrandom::gaussian() + mm::pinning_field_x[cell];
      y_total_external_field_array[cell] = mm::ext_field[1] + sigma_perp*mtrandom::gaussian() + mm::pinning_field_y[cell];
      z_total_external_field_array[cell] = mm::ext_field[2] + sigma_perp*mtrandom::gaussian() + mm::pinning_field_z[cell];

   //   std::cout << pinning_field_y[cell] <<std::endl;
     // optionally add dipole field
      if (dipole::activated){
         x_total_external_field_array[cell] += dipole::cells_field_array_x[cell];
         y_total_external_field_array[cell] += dipole::cells_field_array_y[cell];
         z_total_external_field_array[cell] += dipole::cells_field_array_z[cell];
      }

      // // optionally add bias magnet fields
      if (bias_magnets == true){
         x_total_external_field_array[cell] += bias_field_x[cell];
         y_total_external_field_array[cell] += bias_field_y[cell];
         z_total_external_field_array[cell] += bias_field_z[cell];
      }

      // add dipole field from environment
      if (environment::enabled){
         x_total_external_field_array[cell] += environment::environment_field_x[cell];
         y_total_external_field_array[cell] += environment::environment_field_y[cell];
         z_total_external_field_array[cell] += environment::environment_field_z[cell];
         //std::cout << environment::environment_field_x[cell] << '\t' << environment::environment_field_y[cell] << '\t' << environment::environment_field_z[cell] << std::endl;
      }

      // add field from tracks
      if (sim::track_field_x.size() != 0 ){
         x_total_external_field_array[cell] += sim::track_field_x[cell];
         y_total_external_field_array[cell] += sim::track_field_y[cell];
         z_total_external_field_array[cell] += sim::track_field_z[cell];
         //std::cout << sim::track_field_x[cell] << '\t' << sim::track_field_y[cell] << '\t' << sim::track_field_z[cell] << '\t' << std::endl;
      }

   } // end of loop over all micromagnetic cells

   return;

}
