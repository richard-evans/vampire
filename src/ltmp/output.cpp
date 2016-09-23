//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>

// Vampire headers
#include "ltmp.hpp"
#include "vmpi.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Local variables for writing temperature data
      //-----------------------------------------------------------------------------
      std::ofstream vertical_temperature_file;
      std::ofstream lateral_temperature_file;
      int temperature_profile_output_counter;

      //-----------------------------------------------------------------------------
      // Function declarations
      //-----------------------------------------------------------------------------
      void write_vertical_temperature_data();
      void write_lateral_temperature_data();
      //void write_lateral_vertical_temperature_data(); // To be implemeneted


      //-----------------------------------------------------------------------------
      // Function to output microcell properties
      //-----------------------------------------------------------------------------
      void write_microcell_data(){

         using ltmp::internal::cell_position_array;
         using ltmp::internal::attenuation_array;

         // only output on root process
         if(vmpi::my_rank==0){
            std::ofstream ofile;
            ofile.open("ltmp_cell_coords.cfg");

            for(unsigned int cell=0; cell<cell_position_array.size()/3; ++cell){
               ofile << cell_position_array[3*cell+0] << "\t" << cell_position_array[3*cell+1] << "\t" << cell_position_array[3*cell+2] << "\t";
               ofile << attenuation_array[cell] << std::endl;
            }

            ofile.close();

         }

         return;

      }

      //-----------------------------------------------------------------------------
      // Function to open output file for vertical temperature profile
      //-----------------------------------------------------------------------------
      void open_vertical_temperature_profile_file(){

         temperature_profile_output_counter = 0;
         vertical_temperature_file.open("vertical_temperature_profile.dat");

         return;

      }

      //-----------------------------------------------------------------------------
      // Function to open output file for lateral temperature profile
      //-----------------------------------------------------------------------------
      void open_lateral_temperature_profile_file(){

         temperature_profile_output_counter = 0;
         lateral_temperature_file.open("lateral_temperature_profile.dat");

         return;

      }

      //-----------------------------------------------------------------------------
      // Wrapper function to determine microcell temperature writing function
      //-----------------------------------------------------------------------------
      void write_cell_temperature_data(){

         using ltmp::internal::lateral_discretisation;
         using ltmp::internal::vertical_discretisation;

         if(!lateral_discretisation && vertical_discretisation) ltmp::internal::write_vertical_temperature_data();
         //if(lateral_discretisation && vertical_discretisation) ltmp::internal::write_lateral_vertical_temperature_data;
         if(lateral_discretisation && !vertical_discretisation) ltmp::internal::write_lateral_temperature_data();

         // increment counter
         temperature_profile_output_counter++;

         return;

      }

      //-----------------------------------------------------------------------------
      // Function to write vertical temperature profile to file
      //-----------------------------------------------------------------------------
      void write_vertical_temperature_data(){

         using ltmp::internal::root_temperature_array;

         // only output on root process
         if(vmpi::my_rank==0){
            vertical_temperature_file << temperature_profile_output_counter << "\t";
            for(unsigned int cell=0; cell<root_temperature_array.size()/2; ++cell){
               vertical_temperature_file << root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0] << "\t"; //Te
               vertical_temperature_file << root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] << "\t"; // Tp
            }
            vertical_temperature_file << std::endl;
         }

         return;

      }

      //-----------------------------------------------------------------------------
      // Function to write lateral temperature profile to file
      //-----------------------------------------------------------------------------
      void write_lateral_temperature_data(){

         using ltmp::internal::root_temperature_array;

         // only output on root process
         if(vmpi::my_rank==0){
            lateral_temperature_file << temperature_profile_output_counter << "\t";
            for(unsigned int cell=0; cell<root_temperature_array.size()/2; ++cell){
               lateral_temperature_file << root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0] << "\t"; //Te
               lateral_temperature_file << root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] << "\t"; // Tp
            }
            lateral_temperature_file << std::endl;
         }

         return;

      }

   } // end of internal namespace
} // end of ltmp namespace
