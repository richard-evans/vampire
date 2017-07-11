//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iomanip>
#include <fstream>
#include <sstream>

// Vampire headers
#include "config.hpp"
#include "material.hpp"
#include "sim.hpp"

// config module headers
#include "internal.hpp"

namespace config{
   namespace internal{

      //---------------------------------------------------------------------
      // Function to write meta data for coordinate and atomic data
      //---------------------------------------------------------------------
      void write_coordinate_meta(){

         std::ofstream scmf; // spin coordinate meta file
         scmf.open("atoms-coords.meta");

         // determine file format
         std::string format_string;

         switch(config::internal::format){

            case config::internal::binary:
               format_string = "binary";
               break;

            case config::internal::text:
               format_string = "text";
               break;

         }

         // Get system date
         time_t rawtime = time(NULL);
         struct tm * timeinfo = localtime(&rawtime);

         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Atomistic coordinates configuration file for vampire V5+"<< std::endl;
         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Date: "<< asctime(timeinfo);
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Format: "<< format_string << std::endl;
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Number of atoms: "<< config::internal::total_output_atoms << std::endl;
         scmf << "#--------------------------------------------" << std::endl;
         scmf << "Number of materials: " << mp::num_materials << std::endl;
         for(int mat=0;mat<mp::num_materials;mat++){
            scmf << mat << "\t" << mp::material[mat].mu_s_SI/9.274e-24 << "\t" << mp::material[mat].element << "\t" <<
            mp::material[mat].name << std::endl;
         }
         scmf << "#--------------------------------------------" << std::endl;
         scmf << "Number of coord files: " << config::internal::num_io_groups << std::endl;

         // set simple file name for single file output
         if(config::internal::num_io_groups == 1) scmf << "atoms-coords.data" << std::endl;
         // otherwise set indexed files
         else{
            for(int fid = 0; fid < config::internal::num_io_groups; fid++){
               scmf << "atoms-coords-" << std::setfill('0') << std::setw(6) << fid << ".data" << "\n";
            }
            // flush data to disk
            scmf << std::flush;
         }

         // number of cell files + file list

         // close file
         scmf.close();

      }

      //---------------------------------------------------------------------
      // Function to write meta data for each configuration
      //---------------------------------------------------------------------
      void write_meta(const double simulation_time, // time (seconds)
                      const double temperature, // system temperature (Kelvin)
                      const double applied_field_x, // applied field components (Tesla)
                      const double applied_field_y,
                      const double applied_field_z,
                      const double applied_field_mag,
                      const double magnetization_x, // magnetization components (normalized)
                      const double magnetization_y,
                      const double magnetization_z){

         // determine file name
         std::stringstream filename;
         filename << "spins-";
         filename << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter;
         filename << ".meta";

         // Get system date
         time_t rawtime = time(NULL);
         struct tm * timeinfo = localtime(&rawtime);

         // Declare and open output file
         std::ofstream ofile;
         ofile.open (filename.str());

         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "# Atomistic spin configuration file for vampire v5+"<< "\n";
         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "# Date: "<< asctime(timeinfo);
         ofile << "#------------------------------------------------------"<< "\n";
         ofile << "Time: " << simulation_time << "\n";
         ofile << "Field: " << applied_field_x*applied_field_mag << "\t" << applied_field_y*applied_field_mag << "\t" << applied_field_z*applied_field_mag << "\n";
         ofile << "Temperature: "<< temperature << "\n";
         ofile << "Magnetisation: " << magnetization_x << "\t" << magnetization_y << "\t" << magnetization_z << "\n";
         ofile << "#------------------------------------------------------" << "\n";
         ofile << "Number of spin files: " << config::internal::num_io_groups << "\n"; //vmpi::num_processors-1 << "\n";

         // set simple file name for single file output
         if(config::internal::num_io_groups == 1) ofile << "spins-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".data" << std::endl;
         // otherwise set indexed files
         else{
            for(int fid = 0; fid < config::internal::num_io_groups; fid++){
               ofile << "spins-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << "-" << std::setfill('0') << std::setw(6) << fid << ".data" << "\n";
            }
            // flush data to disk
            ofile << std::flush;
         }

         ofile << "#------------------------------------------------------"<< "\n";

         return;

      }

      //---------------------------------------------------------------------
      // Function to write meta data for coordinate and atomic data
      //---------------------------------------------------------------------
      void write_non_magnetic_meta(const uint64_t num_data){

         std::ofstream scmf; // spin coordinate meta file
         scmf.open("non-magnetic-atoms.meta");

         // determine file format
         std::string format_string;

         switch(config::internal::format){

            case config::internal::binary:
               format_string = "binary";
               break;

            case config::internal::text:
               format_string = "text";
               break;

         }

         // Get system date
         time_t rawtime = time(NULL);
         struct tm * timeinfo = localtime(&rawtime);

         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Atomistic coordinates configuration file for vampire V5+"<< std::endl;
         scmf << "#----------------------------------------------------------"<< std::endl;
         scmf << "# Date: "<< asctime(timeinfo);
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Format: "<< format_string << std::endl;
         scmf << "#--------------------------------------------"<< std::endl;
         scmf << "Number of atoms: "<< num_data << std::endl;
         scmf << "#--------------------------------------------" << std::endl;
         scmf << "Number of files: " << config::internal::num_io_groups << std::endl;

         // set simple file name for single file output
         if(config::internal::num_io_groups == 1) scmf << "non-magnetic-atoms.data" << std::endl;
         // otherwise set indexed files
         else{
            for(int fid = 0; fid < config::internal::num_io_groups; fid++){
               scmf << "non-magnetic-atoms-" << std::setfill('0') << std::setw(6) << fid << ".data" << "\n";
            }
            // flush data to disk
            scmf << std::flush;
         }

         // number of cell files + file list

         // close file
         scmf.close();

      }

   } // end of namespace internal
} // end of namespace config
