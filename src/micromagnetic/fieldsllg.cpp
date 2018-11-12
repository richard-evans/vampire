//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "environment.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// micromagnetic module headers
#include "micromagnetic.hpp"
#include "internal.hpp"

namespace micromagnetic{

namespace internal{

// stack variables for file
std::vector<double> spin_field(3,0.0);
const double twokB = 2.0*1.3806503e-23; // Boltzmann's constant

//------------------------------------------------------------------------------
// Function to calculate field for LLG micromagnetics
//------------------------------------------------------------------------------
std::vector<double> calculate_llg_fields(std::vector <double>& m,
                                         double temperature,
                                         int num_cells,
                                         int cell,
                                         std::vector<double>& x_array,
                                         std::vector<double>& y_array,
                                         std::vector<double>& z_array){


   //the temperature is usually used as a reduced temperature.
   const double reduced_temperature = temperature/Tc[cell];

   // calculate equilibrium magnetization m_e
   if(reduced_temperature <= 1.0){
      m_e[cell] = pow((1.0 - reduced_temperature),0.34); // standard fit to Ms(T) for classical model
   }
   else{
      // assume small minimum value to avoid numerical issues above Tc for LLG
      m_e[cell] = 0.01;
   }


   //------------------------------------------------------------------------------
   // Calculate the exchange fields as me^1.66 *A*(xi-xj)/m_e^2
   //------------------------------------------------------------------------------

   // temporary variable to hold exchange field
   double exchange_field[3]={0.0,0.0,0.0};

   // temporary constants for cell i
   const double mi = m_e[cell];
   const double mix = mi * m[0];
   const double miy = mi * m[1];
   const double miz = mi * m[2];
   //const double m_e_squared = mi*mi;

   // If there is more than one micromagnetic cell then calculate exchange, otherwise skip
   if (num_cells > 1){

      const int start = macro_neighbour_list_start_index[cell]; // start of neighbour list
      const int end   = macro_neighbour_list_end_index[cell];   // end of neighbour list

      // loop over all neighbours
      for(int j = start; j<end; j++){

         // calculate reduced exchange constant factor
         const int cellj = macro_neighbour_list_array[j]; // get neighbour cell ID
         const double mj = m_e[cellj]; // get reduced magnetization of neighbour cell
         double Ac = A[j]*pow(mj,1.66);

         // check for SAF interaction between cells
         int mat  = cell_material_array[cell];
         int matj = cell_material_array[cellj];
         if (mp::material[mat].enable_SAF && mp::material[matj].enable_SAF){
            if (mat != matj){
               double Area = cells::macro_cell_size[0]*cells::macro_cell_size[1];
               Ac = -pow(mj,1.66)*Area*mp::material[mat].SAF[matj]/ms[cell];
               //if (mm_correction == true) Ac = 2*Ac/cells::macro_cell_size[2]; // what does this do?
            }
         }

         // add exchange field to sum
         exchange_field[0] -= Ac * (x_array[cellj]*mj - mix);
         exchange_field[1] -= Ac * (y_array[cellj]*mj - miy);
         exchange_field[2] -= Ac * (z_array[cellj]*mj - miz);

      }

   }

   // calculate sigma value (added fixed alpha and m_e in denominator?)
   const double sigma = sqrt( twokB * temperature * alpha[cell] / ( mp::dt * ms[cell] * mi ) );

   const double anis_field[3]    = {one_o_chi_perp[cell]*mix, one_o_chi_perp[cell]*miy, 0.0};
   const double thermal_field[3] = {sigma * mtrandom::gaussian(), sigma * mtrandom::gaussian(), sigma * mtrandom::gaussian()};

   spin_field[0] = ext_field[0] + exchange_field[0] + anis_field[0] + thermal_field[0] + pinning_field_x[cell];
   spin_field[1] = ext_field[1] + exchange_field[1] + anis_field[1] + thermal_field[1] + pinning_field_y[cell];
   spin_field[2] = ext_field[2] + exchange_field[2] + anis_field[2] + thermal_field[2] + pinning_field_z[cell];

   // Add dipole field if enabled
   if (dipole::activated){
      spin_field[0] += dipole::cells_field_array_x[cell];
      spin_field[1] += dipole::cells_field_array_y[cell];
      spin_field[2] += dipole::cells_field_array_z[cell];
   }

   // Add FMR field if enabled
   if (sim::enable_fmr){
      spin_field[0] += fmr_H[0];
      spin_field[1] += fmr_H[1];
      spin_field[2] += fmr_H[2];
   }

   // Add field from environment if enabled
   if (environment::enabled){
      spin_field[0] += environment::environment_field_x[cell];
      spin_field[1] += environment::environment_field_y[cell];
      spin_field[2] += environment::environment_field_z[cell];
   }

   // Add field from tracks if enabled
   if (sim::track_field_x.size() != 0 ){
      spin_field[0] += sim::track_field_x[cell];
      spin_field[1] += sim::track_field_y[cell];
      spin_field[2] += sim::track_field_z[cell];
   }

   return spin_field;

}

} // end of internal namespace

} // end of micromagnetic namespace
