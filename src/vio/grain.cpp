//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "vio.hpp"
#include "grains.hpp"
#include "sim.hpp"
#include "stats.hpp"

// vio module headers
#include "internal.hpp"
#include "stats.hpp"

namespace vout{

void write_grain_file(){

   // do nothing if no grain items specified
   if(vout::grain::output_list.size() == 0) return;

   // check it is time to output a new data point
   if(sim::time % vout::grain::output_rate == 0){

      // disable headers for variables
      bool header = false;

      // Output data to zgrain
      if(vmpi::my_rank==0){

         // check for open ofstream
         if( !zgrain.is_open() ){
            // check for checkpoint continue and append data
            if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) zgrain.open("grain.txt",std::ofstream::app);
            // otherwise overwrite file
            else zgrain.open("grain.txt",std::ofstream::trunc);
         }

         // loop over all data in the output data list
         for(unsigned int item=0; item < vout::grain::output_list.size(); item++){

            switch(vout::grain::output_list[item]){
               //------------------------------------------
               case grain::time_steps:
                  vout::time(zgrain,header);
                  break;
               //------------------------------------------
               case grain::real_time:
                  vout::real_time(zgrain,header);
                  break;
               //------------------------------------------
               case grain::temperature:
                  vout::temperature(zgrain,header);
                  break;
               //------------------------------------------
               case grain::electron_temperature:
                  vout::temperature(zgrain,header);
                  break;
               //------------------------------------------
               case grain::phonon_temperature:
                  vout::phonon_temperature(zgrain,header);
                  break;
               //------------------------------------------
               case grain::applied_field:
                  vout::Happ(zgrain,header);
                  break;
               //------------------------------------------
               case grain::applied_field_unit_vector:
                  vout::Hvec(zgrain,header);
                  break;
               //------------------------------------------
               case grain::constraint_phi:
                  zgrain << generic_output_double("con_phi",sim::constraint_phi,header);
                  break;
               //------------------------------------------
               case grain::constraint_theta:
                  zgrain << generic_output_double("con_theta",sim::constraint_theta,header);
                  break;
               //------------------------------------------
               case grain::magnetisation:
                  // inline function to output grain data
                  zgrain << stats::grain_magnetization.output_normalized_magnetization(header);
                  break;
               //------------------------------------------
               case grain::mean_magnetisation_length:
                  // inline function to output grain data
                  zgrain << stats::grain_magnetization.output_normalized_mean_magnetization_length(header);
                  break;
               //------------------------------------------
               case grain::material_magnetisation:
                  // inline function to output grain data
                  zgrain << stats::material_grain_magnetization.output_normalized_magnetization(header);
                  break;
               //------------------------------------------
               case grain::material_height_magnetisation:
                  // inline function to output grain data
                  zgrain << stats::material_grain_height_magnetization.output_normalized_magnetization(header);
                  break;
               //------------------------------------------
               case grain::mean_torque:
                  // inline function to output grain data
                  zgrain << stats::grain_torque.output_mean_torque(header);
                  break;
               //------------------------------------------
               case grain::mean_susceptibility:
                  // inline function to output grain data
                  zgrain << stats::grain_susceptibility.output_mean_susceptibility(sim::temperature, header);
                  break;
               //------------------------------------------
               case grain::mean_specific_heat:
                  // inline function to output grain data
                  zgrain << stats::grain_specific_heat.output_mean_specific_heat(sim::temperature, header);
                  break;
               //------------------------------------------

            } // end of case statement

         } // end of output list loop

         // Carriage return
         zgrain << std::endl;

      } // end of rank zero check

   } // end of output check

   return;

} // end of function

} // end of namespace vout
