//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
bool read_spin_metadata(unsigned int file_id);
void read_spin_data();

//------------------------------------------------------------------------------
// Wrapper function to read coordinate metafile to initialise data structures
// and process coordinate data
//------------------------------------------------------------------------------
void process_spins(){

   unsigned int min_file_id = vdc::vdc_start_file_id;
   unsigned int max_file_id = vdc::vdc_final_file_id;

   if(vdc::cells) vdc::initialise_cells();

   if(vdc::ssc) vdc::initialise_ssc();

   if(vdc::povray || vdc::povcells ) vdc::initialise_povray();

   // output povray file
   if(vdc::povray) output_povray_file();
   if(vdc::povcells) output_povray_cells_file();

   unsigned int last_file_id = max_file_id;

   // loop over all spin files
   for(unsigned int file_id = min_file_id; file_id < max_file_id; file_id++){

      // read meta data
      bool success = vdc::read_spin_metadata(file_id);

      // if no success then break out of for loop
      if(!success) break;

      // read coordinate data
      vdc::read_spin_data();

      // output cells raw data
      if(vdc::cells) vdc::output_cell_file(file_id);

      // output povray files
      if(vdc::povray) output_inc_file(file_id);
      if(vdc::povcells) output_cells_inc_file(file_id);

      // output vtk file
      if(vdc::vtk) output_vtk_file(file_id);

      // output plain text file
      if(vdc::txt) output_txt_file(file_id);

      // compute spin-spin correlation
      if(vdc::ssc) output_ssc_file(file_id);

      last_file_id = file_id;

   }

   // set global start and end file id
   vdc::start_file_id = min_file_id;
   vdc::final_file_id = last_file_id;

   // output average ssc
   if(vdc::ssc) output_average_ssc_file();

   return;

}

//------------------------------------------------------------------------------
// Function to read coordinate metafile
//------------------------------------------------------------------------------
//
// Example metafile format:
//
//       #------------------------------------------------------
//       # Atomistic spin configuration file for vampire v5+
//       #------------------------------------------------------
//       # Date: Fri Apr  7 21:42:33 2017
//       #------------------------------------------------------
//       Time: 1e-14
//       Field: 0	0	1
//       Temperature: 0
//       Magnetisation: 0.707068	5.51664e-05	0.707146
//       #------------------------------------------------------
//       Number of spin files: 1
//       spins-00000000.data
//       #------------------------------------------------------
//
//------------------------------------------------------------------------------
bool read_spin_metadata(unsigned int file_id){

   // determine file name
   std::stringstream filename;
   filename << "spins-";
   filename << std::setfill('0') << std::setw(8) << file_id;
   filename << ".meta";

   // open spins metadata file
   std::ifstream smfile;
   smfile.open(filename.str());

   // check for open file, if not open then end program, end of snapshots
   if(!smfile.is_open()){
         //std::cerr << "Error! Spins metadata file spins-" << std::setfill('0') << std::setw(8)
         //<< file_id << ".meta cannot be opened. Exiting" << std::endl;
         //exit(1);
      return false;
   }

   // Metafile found - inform the user and process data
   if(vdc::verbose) std::cout << "--------------------------------------------------------------------" << std::endl;
   std::cout << "Processing snapshot " << std::setfill('0') << std::setw(8) << file_id << std::endl;
   if(vdc::verbose) std::cout << "   Reading spin meta-data file " << filename.str() << std::endl;

   std::string line; // line string variable

   // read in file header (not useful - need to read in variables)
   for(int i=0; i<10; i++) getline(smfile, line);

   // get number of subsidiary files
   getline(smfile, line);
   line.erase (line.begin(), line.begin()+22);
   unsigned int num_spin_files=atoi(line.c_str());

   if(vdc::verbose) std::cout << "   Number of data files: " << num_spin_files << std::endl;

   vdc::spin_filenames.resize(0);

   for(unsigned int file = 0; file < num_spin_files; file++){
      getline(smfile, line);
      line.erase(remove(line.begin(), line.end(), '\t'), line.end());
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      line.erase(remove(line.begin(), line.end(), '\r'), line.end());
      vdc::spin_filenames.push_back(line);
      if(vdc::verbose) std::cout << "      " << line << std::endl;
   }

   return true;

}

//------------------------------------------------------------------------------
// Function to read in coordinate data from subsidiary files
//------------------------------------------------------------------------------
void read_spin_data(){

   if(vdc::verbose) std::cout << "   Reading spin data... " << std::flush;

   // resize arrays
   if(vdc::spins.size() != 3*vdc::num_atoms) vdc::spins.resize(3*vdc::num_atoms);

   // index counter
   uint64_t atom_id = 0;

   // loop over all files
   for(unsigned int f = 0; f < vdc::spin_filenames.size(); f++){

      switch (vdc::format){

         case vdc::binary:{
            uint64_t num_atoms_in_file = 0;
            // open file in binary mode
            std::ifstream ifile;
            ifile.open(spin_filenames[f].c_str(), std::ios::binary); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << std::endl << "   Error! Spin data file \"" << spin_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
               exit(1);
            }
            // read number of atoms
            ifile.read( (char*)&num_atoms_in_file,sizeof(uint64_t) );
            // read spin data
            ifile.read((char*)&vdc::spins[atom_id*3], sizeof(double)*num_atoms_in_file*3);
            // increment counter
            atom_id += num_atoms_in_file;
            ifile.close();
            break;
         }

         case vdc::text:{
            // open file
            std::ifstream ifile;
            ifile.open(spin_filenames[f].c_str()); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << std::endl << "   Error! Spin data file \"" << spin_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
               exit(1);
            }

            uint64_t num_atoms_in_file = 0;
            std::string line;
            getline(ifile, line);
            {
               std::istringstream ss(line);
               ss >> num_atoms_in_file; // interpret as uint64_t
            }
            double x,y,z;
            // loop over all atoms in file and load as x,y,z sets
            for(uint64_t idx = 0; idx < num_atoms_in_file; idx++){
               getline(ifile, line);
               std::istringstream ss(line);
               ss >> x >> y >> z;
               vdc::spins[3*atom_id + 0] = x;
               vdc::spins[3*atom_id + 1] = y;
               vdc::spins[3*atom_id + 2] = z;
               // increment atom counter
               atom_id += 1;
            }
            ifile.close();
            break;
         }

      }

   }

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
