//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2015. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iomanip>
#include <sstream>

// Vampire headers
#include "config.hpp"

// config headers
#include "internal.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions to output configuration files
//--------------------------------------------------------------------------------
namespace config{

   //-----------------------------------------------------------------------------
   //
   //
   //
   //
   //
   //
   //
   //
   //
   //-----------------------------------------------------------------------------
   void output(const std::vector<double>& spins_x, // spin magnetization vectors (unit)
               const std::vector<double>& spins_y,
               const std::vector<double>& spins_z,
               const std::vector<double>& cells_x, // cell magnetization vectors (Bohr magnetons)
               const std::vector<double>& cells_y,
               const std::vector<double>& cells_z,
               const double simulation_time, // time (seconds)
               const double temperature, // system temperature (Kelvin)
               const double applied_field_x, // applied field components (Tesla)
               const double applied_field_y,
               const double applied_field_z,
               const double magnetization_x, // magnetization components (normalized)
               const double magnetization_y,
               const double magnetization_z){

      // Determine if configuration should be output
      bool output_config_now = (config::internal::output_rate_counter % config::internal::output_rate) == 0;
      
      if(output_config_now){

         
         // meta data output
         if(config::internal::output_meta){
            /*config::internal::write_meta(simulation_time, temperature, 
                                         applied_field_x, applied_field_y, applied_field_z,
                                         magnetization_x, magnetization_y, magnetization_z);*/
         }
         
         // atoms output
         if(config::internal::output_atoms){
            
            // determine filename
            std::stringstream filename;
            filename << "atoms-";
            filename << std::setfill('0') << std::setw(8) << config::internal::output_file_counter;
            filename << ".cfg";

            // copy to buffer
            config::internal::copy_data_to_buffer(spins_x, spins_y, spins_z, config::internal::local_output_atom_list, config::internal::output_spin_buffer); 
            
            // write data to disk
            config::internal::write_data(filename.str(),config::internal::output_spin_buffer);
            
         }
         
         // cells output
         //if(config::internal::output_macrocells) write_macrocells(cells_x, cells_y, cells_z);

      }

      // increment rate counter
      config::internal::output_rate_counter++;
         
      return;

   }

} // end of namespace config
