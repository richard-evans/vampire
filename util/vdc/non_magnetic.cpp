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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

//------------------------------------------------------------------------------
// Function to read coordinate metafile
//------------------------------------------------------------------------------
//
// Example metafile format:
//
//    #----------------------------------------------------------
//    # Atomistic coordinates configuration file for vampire V5+
//    #----------------------------------------------------------
//    # Date: Tue Jun 27 23:01:08 2017
//    #--------------------------------------------
//    Format: text
//    #--------------------------------------------
//    Number of atoms: 3564
//    #--------------------------------------------
//    Number of files: 1
//    non-magnetic-atoms.data
//
//------------------------------------------------------------------------------
void read_nm_metadata(){

   if(vdc::verbose) std::cout << "Reading non-magnetic meta-data" << std::endl;

   // open coordinate metadata file
   std::ifstream cmfile;
   cmfile.open("non-magnetic-atoms.meta");

   // check for open file
   if(!cmfile.is_open()){
      return;
   }

   std::string line; // line string variable

   // read in file header (not useful)
   for(int i=0; i<5; i++) getline(cmfile, line);

   // get data format
   getline(cmfile, line);
   line.erase (line.begin(), line.begin()+8);
   std::string data_format_str = line;

   line.erase(remove(line.begin(), line.end(), '\t'), line.end());
   line.erase(remove(line.begin(), line.end(), ' '), line.end());
   line.erase(remove(line.begin(), line.end(), '\r'), line.end());

   // Tell the user about the data format if verbosity is required
   if(vdc::verbose) std::cout << "   Data format: " << data_format_str << std::endl;

   std::string test = "text";
   if(data_format_str == test){
     vdc::format = vdc::text;
     if(vdc::verbose) std::cout << "   Setting data format to text mode" << std::endl;
   }
   test = "binary";
   if(data_format_str == test){
     vdc::format = vdc::binary;
     if(vdc::verbose) std::cout << "   Setting data format to binary mode" << std::endl;
   }
   /*else{
      std::cerr << "Unknown data format \"" << data_format_str << "\". Exiting" << std::endl;
      exit(1);
   }*/

   // get comment line
   getline(cmfile, line);

   // get number of atoms
   getline(cmfile, line);
   line.erase (line.begin(), line.begin()+17);
   {
      std::istringstream ss(line);
      ss >> vdc::num_nm_atoms; // interpret as uint64_t
   }

   // Tell the user about the total number of atoms if verbosity is required
   if(vdc::verbose) std::cout << "   Number of non-magnetic atoms: " << num_nm_atoms << std::endl;

   // get comment line
   getline(cmfile, line);

   // get number of subsidiary files
   unsigned int num_coord_files;
   getline(cmfile, line);
   line.erase (line.begin(), line.begin()+17);
   num_coord_files=atoi(line.c_str());

   if(vdc::verbose) std::cout << "   Number of files: " << num_coord_files << std::endl;

   for(unsigned int file = 0; file < num_coord_files; file++){
      getline(cmfile, line);
      line.erase(remove(line.begin(), line.end(), '\t'), line.end());
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      line.erase(remove(line.begin(), line.end(), '\r'), line.end());
      vdc::nm_filenames.push_back(line);
      if(vdc::verbose) std::cout << "      " << line << std::endl;
   }

   return;

}

//------------------------------------------------------------------------------
// Function to read in coordinate data from subsidiary files
//------------------------------------------------------------------------------
void read_nm_data(){

   if(vdc::verbose) std::cout << "Reading non-magnetic data... " << std::flush;

   // resize arrays
   vdc::nm_coordinates.resize(3*vdc::num_nm_atoms);
   vdc::nm_category.resize(vdc::num_nm_atoms);
   vdc::nm_type.resize(vdc::num_nm_atoms);

   // index counter
   uint64_t atom_id = 0;

   // loop over all files
   for(unsigned int f = 0; f < vdc::nm_filenames.size(); f++){

      switch (vdc::format){

         case vdc::binary:{
            uint64_t num_atoms_in_file = 0;
            // open file in binary mode
            std::ifstream ifile;
            ifile.open(nm_filenames[f].c_str(), std::ios::binary); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << "Error! non-magnetic data file \"" << nm_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
               std::cerr << "Error code: " << std::strerror(errno) << std::endl;
               exit(1);
            }
            // read number of atoms
            ifile.read( (char*)&num_atoms_in_file,sizeof(uint64_t) );
            // read type array
            ifile.read((char*)&vdc::nm_type[atom_id], sizeof(int)*num_atoms_in_file);
            // read category array
            ifile.read((char*)&vdc::nm_category[atom_id], sizeof(int)*num_atoms_in_file);
            ifile.read((char*)&vdc::nm_coordinates[3*atom_id], sizeof(double)*num_atoms_in_file*3);
            // increment counter
            atom_id += num_atoms_in_file;
            ifile.close();

            break;
         }

         case vdc::text:{
            // open file
            std::ifstream ifile;
            ifile.open(nm_filenames[f].c_str()); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << "Error! non-magnetic data file \"" << nm_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
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
            int type_id, category_id;
            // loop over all atoms in file and load as x,y,z sets
            for(uint64_t idx = 0; idx < num_atoms_in_file; idx++){
               getline(ifile, line);
               std::istringstream ss(line);
               ss >> type_id >> category_id >> x >> y >> z;
               vdc::nm_type[atom_id] = type_id;
               vdc::nm_category[atom_id] = category_id;
               vdc::nm_coordinates[3*atom_id + 0] = x;
               vdc::nm_coordinates[3*atom_id + 1] = y;
               vdc::nm_coordinates[3*atom_id + 2] = z;
               // increment atom counter
               atom_id += 1;
            }
            ifile.close();
            break;
         }

      }

   }

   // create vector list of non magnetic atom indices
   vdc::nm_atoms_list.resize(vdc::num_nm_atoms);
   for(unsigned int atom = 0; atom < vdc::num_nm_atoms; atom++){
      vdc::nm_atoms_list.push_back(atom);
   }

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   // check that number of materials is sufficient
   int num_materials = 0;
   for(int mat : vdc::nm_type){
      if(mat > num_materials-1) num_materials = mat+1;
   }
   if(static_cast<unsigned int>(num_materials) > materials.size()) materials.resize(num_materials);

   return;

}

}
