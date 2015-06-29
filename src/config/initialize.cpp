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
#include "material.hpp"
#include "vio.hpp"

// config headers
#include "internal.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions to output configuration files
//--------------------------------------------------------------------------------
namespace config{

   // Forward declaration of functions
   
   
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
   void initialize(const int num_atoms, // number of local atoms
                   const double system_dimensions_x, // system size
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<double>& spins_cx, // spin coordinates (Angstroms)
                   const std::vector<double>& spins_cy,
                   const std::vector<double>& spins_cz,
                   const std::vector<double>& cells_cx, // cell coordinates (Angstroms)
                   const std::vector<double>& cells_cy,
                   const std::vector<double>& cells_cz,
                   const std::vector<int>& spins_mat)  // spin material ID
   {

      // Check if config output is needed
      if(config::internal::output_meta || config::internal::output_coords){
         
      zlog << zTs() << "Initializing configuration file output" << std::endl;

      // determine number of atoms to output on io process
      // mpi_all_reduce(num_atoms) MPI_COMM_NODE
      
      // allocate temporary storage on vmpi::my_io
      
      //--------------------------------------------------------------------
      // Determine atoms to output
      //--------------------------------------------------------------------
      
      // resize atom list to zero
      config::internal::local_output_atom_list.resize(0);
      
      // get output bounds
      const double minB[3]={config::internal::atoms_output_min[0]*system_dimensions_x,
                            config::internal::atoms_output_min[1]*system_dimensions_y,
                            config::internal::atoms_output_min[2]*system_dimensions_z};

      const double maxB[3]={config::internal::atoms_output_max[0]*system_dimensions_x,
                            config::internal::atoms_output_max[1]*system_dimensions_y,
                            config::internal::atoms_output_max[2]*system_dimensions_z};

      // loop over all local atoms and record output list
      for(int atom=0;atom<num_atoms;atom++){

         const double cc[3] = {spins_cx[atom],spins_cy[atom],spins_cz[atom]};

         // check atom within output bounds
         if((cc[0] >= minB[0]) && (cc[0]<=maxB[0])){
            if((cc[1] >= minB[1]) && (cc[1]<=maxB[1])){
               if((cc[2] >= minB[2]) && (cc[2]<=maxB[2])){
                  config::internal::local_output_atom_list.push_back(atom);
               }
            }
         }
      }
      
      config::internal::total_output_atoms = config::internal::local_output_atom_list.size();
      
      // resize output buffer to 3*num_output_atoms
      config::internal::output_spin_buffer.resize(3*config::internal::total_output_atoms);
      // cell buffer...
      
      //-------------------------------------------------------
      // Output spin coordinate meta data
      //-------------------------------------------------------
      if(config::internal::output_atoms || config::internal::output_coords){
      
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

            //

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
            scmf << "Number of coord files: " << 1 << std::endl;
            scmf << "atoms-coords.cfg" << std::endl;

            // number of cell files + file list
            
            // close file
            scmf.close();
            
            // output coordinate data
            
            // determine filename
            std::stringstream filename;
            filename << "atoms-coords.cfg";

            // copy to buffer
            config::internal::copy_data_to_buffer(spins_cx, spins_cy, spins_cz, config::internal::local_output_atom_list, config::internal::output_spin_buffer); 
            
            // write buffer to disk
            config::internal::write_data(filename.str(),config::internal::output_spin_buffer);      
            
         #endif
         
      }
      
      }
      
      return;

   }

   
   
} // end of namespace config
