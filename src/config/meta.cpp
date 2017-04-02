//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
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

         #ifdef MPICF

         #else
            std::ofstream scmf; // spin coordinate meta file
            scmf.open("atoms-coords.meta");

            // determine file format
            std::string format_string;

            switch(config::internal::output_data_format){

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
            scmf << "Number of atoms: "<< config::total_output_atoms << std::endl;
            scmf << "#--------------------------------------------" << std::endl;
            scmf << "Number of materials: " << mp::num_materials << std::endl;
            for(int mat=0;mat<mp::num_materials;mat++){
               scmf << mat << "\t" << mp::material[mat].mu_s_SI/9.274e-24 << "\t" << mp::material[mat].element << "\t" <<
               mp::material[mat].name << std::endl;
            }
            scmf << "#--------------------------------------------" << std::endl;
            scmf << "Number of coord files: " << 1 << std::endl;
            scmf << "atoms-coords.cfg" << std::endl;

            // number of cell files + file list

            // close file
            scmf.close();

         #endif

      }

      //---------------------------------------------------------------------
      // Function to write meta data for each configuration
      //---------------------------------------------------------------------
      void write_meta(const double simulation_time, // time (seconds)
                      const double temperature, // system temperature (Kelvin)
                      const double applied_field_x, // applied field components (Tesla)
                      const double applied_field_y,
                      const double applied_field_z,
                      const double magnetization_x, // magnetization components (normalized)
                      const double magnetization_y,
                      const double magnetization_z,
                      const int    num_files){

         // determine file name
         std::stringstream filename;
         filename << "atoms-";
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
         ofile << "Field: " << applied_field_x << "\t" << applied_field_y << "\t" << applied_field_z << "\n";
         ofile << "Temperature: "<< temperature << "\n";
         ofile << "Magnetisation: " << magnetization_x << "\t" << magnetization_y << "\t" << magnetization_z << "\n";
         ofile << "#------------------------------------------------------" << "\n";
         ofile << "Number of spin files: " << num_files << "\n"; //vmpi::num_processors-1 << "\n";
         for(int p=0;p<num_files;p++){
            std::stringstream cfg_sstr;


            cfg_sstr << "atom-spins-" ;
            switch (config::internal::format)
            {
            case config::internal::binary:
               cfg_sstr << "binary-";
               break;
            case config::internal::text:
               cfg_sstr << "text-";
               break;
            }
            cfg_sstr << std::setfill('0') << std::setw(5) << p << "-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".cfg";
            ofile << cfg_sstr.str() << "\n";
         }

         ofile << "#------------------------------------------------------"<< "\n";

         return;

      }

   } // end of namespace internal
} // end of namespace config
