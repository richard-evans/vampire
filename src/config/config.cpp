//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo, Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp"
#include "config.hpp"
#include "gpu.hpp"
#include "program.hpp"
#include "sim.hpp"

// config module headers
#include "internal.hpp"

namespace config
{

//------------------------------------------------------------------------------
// Function to output atomic and cell coordinates to disk
//------------------------------------------------------------------------------
void output(){ // should include variables for data to be outputted, eg spins, cells etc

   // ALWAYS enable output of atomic configurations at end of simulation for HAMR simulation, even if output wasn't requested
   if (program::program == 7){
      // overwrite these variables
      config::internal::output_atoms_config = true;
      config::internal::output_atoms_config_end = true;
   }

   // check for data output enabled, if not no nothing
   if(config::internal::output_atoms_config == false && config::internal::output_cells_config == false) return;

   // check that config module has been initialised
   if(!config::internal::initialised) config::internal::initialize();

   //------------------------------------------------------------------------------------------
   // Calculate field ranges for output in limited applied ranges during hysteresis
   //------------------------------------------------------------------------------------------

   // Max and Min field values for descending branch
   double minField_1;
   double maxField_1;
   // Max and Min field values for ascending branch
   double minField_2;
   double maxField_2;

   // check that minField_1<maxField_1
   if (config::internal::field_output_min_1 < config::internal::field_output_max_1)
   {
      minField_1 = config::internal::field_output_min_1;
      maxField_1 = config::internal::field_output_max_1;
   }
   else
   {
      minField_1 = config::internal::field_output_max_1;
      maxField_1 = config::internal::field_output_min_1;
   }
   // check that maxField_2>minField_2
   if (config::internal::field_output_max_2 >= config::internal::field_output_min_2)
   {
      minField_2 = config::internal::field_output_min_2;
      maxField_2 = config::internal::field_output_max_2;
   }
   else
   {
      minField_2 = config::internal::field_output_max_2;
      maxField_2 = config::internal::field_output_min_2;
   }

   //------------------------------------------------------
   // atoms output if enabled and the time is right
   //------------------------------------------------------
   if (config::internal::output_atoms_config==true){
      // Print at the end of simulation
      if (config::internal::output_atoms_config_end==true && (sim::time ==  sim::total_time+sim::equilibration_time))
      {
         // If using GPU acceleration then synchonise spins from device
         gpu::config::synchronise();

         //Always output coordinates the first time (re-)started, otherwise the spins coordinates won't be printed
         if (config::internal::output_rate_counter_coords == 0)
         {
             config::internal::atoms_coords();
             if(atoms::num_non_magnetic_atoms > 0) config::internal::atoms_non_magnetic();
          }
         config::internal::atoms(); // call function to output spins coords
         config::internal::output_rate_counter_coords++; //update the counter
      }
      // Otherwise print during simulation
      else if (config::internal::output_atoms_config_continuous && (sim::output_rate_counter % config::internal::output_atoms_config_rate == 0))
      {

         // If using GPU acceleration then synchonise spins from device
         gpu::config::synchronise();

         // for all programs except hysteresis(=2), static-hysteresis(=3) and partial-hysteresis(=12)
         if ((program::program != 2) && (program::program != 3) && (program::program != 12))
         {

            //Always output coordinates the first time (re-)started, otherwise the spins coordinates won't be printed
            if (config::internal::output_rate_counter_coords == 0)
            {
                config::internal::atoms_coords();
                if(atoms::num_non_magnetic_atoms > 0) config::internal::atoms_non_magnetic();
             }
            config::internal::atoms(); // call function to output spins coords
            config::internal::output_rate_counter_coords++; //update the counter
         }
         // for hysteresis programs
         else if ((program::program == 2) || (program::program ==3) || (program::program ==12))
         {
            // output config only in range [minField_1;maxField_1] for descending branch
            if (sim::parity < 0)
            {
               if((sim::H_applied >= minField_1) && (sim::H_applied <= maxField_1))
               {
                  if(config::internal::output_rate_counter_coords == 0)
                  {
                     config::internal::atoms_coords();
                     if(atoms::num_non_magnetic_atoms > 0) config::internal::atoms_non_magnetic();
                  }
                  config::internal::atoms();
                  config::internal::output_rate_counter_coords++;
               }
            }
            // output config only in range [minField_2;maxField_2] for ascending branch
            else if (sim::parity > 0)
            {
               if((sim::H_applied >= minField_2) && (sim::H_applied <= maxField_2)){
                  if (config::internal::output_rate_counter_coords == 0)
                  {
                     config::internal::atoms_coords();
                     if(atoms::num_non_magnetic_atoms > 0) config::internal::atoms_non_magnetic();
                  }
                  config::internal::atoms();
                  config::internal::output_rate_counter_coords++;
               }
            }
         }
      }
   } // if for atoms output

   //------------------------------------------------------
   // cells output if enabled and the time is right
   //------------------------------------------------------
   if (config::internal::output_cells_config){
      // Print only at the end of simulation
      if(config::internal::output_cells_config_end && (sim::time ==  sim::total_time+sim::equilibration_time))
      {
         config::internal::legacy_cells_coords();
         config::internal::legacy_cells();
      }
      // Otherwise print cells configurations during simulation
      else if (config::internal::output_cells_config_continuous  && (sim::output_rate_counter % config::internal::output_cells_config_rate == 0))
      {
         // for all programs except hysteresis(=2), static-hysteresis(=3) and partial-hysteresis(=12)
         if ((program::program != 2) && (program::program != 3) && (program::program != 12))
         {
            if (sim::output_cells_file_counter == 0) config::internal::legacy_cells_coords();
            config::internal::legacy_cells();
         }
         // for hysteresis programs
         else if ((program::program == 2) || (program::program ==3) || (program::program ==12))
         {
            // output config only in range [minField_1;maxField_1] for decreasing field
            if (sim::parity < 0)
            {
               if((sim::H_applied >= minField_1) && (sim::H_applied <= maxField_1))
               {
                  if (sim::output_cells_file_counter == 0) config::internal::legacy_cells_coords();
                  config::internal::legacy_cells();
               }
            }
            // output config only in range [minField_2;maxField_2] for increasing field
            else if (sim::parity > 0)
            {
               if((sim::H_applied >= minField_2) && (sim::H_applied <= maxField_2))
               {
                  if (sim::output_cells_file_counter == 0) config::internal::legacy_cells_coords();
                  config::internal::legacy_cells();
               }
            }
         }
      }
   } //end output of cells

   // increment rate counter
   sim::output_rate_counter++;
}

}
