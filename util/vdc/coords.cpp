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

// forward function declarations
void read_coord_metadata();
void read_coord_data();

//------------------------------------------------------------------------------
// Wrapper function to read coordinate metafile to initialise data structures
// and process coordinate data
//------------------------------------------------------------------------------
void process_coordinates(){

   // read meta data
   vdc::read_coord_metadata();

   // read coordinate data
   vdc::read_coord_data();

   // load non-magnetic data files
   read_nm_metadata();
   read_nm_data();

   // output xyz file
   if(vdc::xyz) output_xyz_file();

   return;

}

//------------------------------------------------------------------------------
// Function to read coordinate metafile
//------------------------------------------------------------------------------
//
// Example metafile format:
//
//    #----------------------------------------------------------
//    # Atomistic coordinates configuration file for vampire V5+
//    #----------------------------------------------------------
//    # Date: Thu Apr  6 23:41:37 2017
//    #--------------------------------------------
//    Format: text
//    #--------------------------------------------
//    Number of atoms: 12167
//    #--------------------------------------------
//    Number of materials: 1
//    0	1.72	Ag 	Co
//    #--------------------------------------------
//    Number of coord files: 1
//    atoms-coords.data
//
//
//------------------------------------------------------------------------------
void read_coord_metadata(){

   if(vdc::verbose) std::cout << "Reading meta-data" << std::endl;

   // open coordinate metadata file
   std::ifstream cmfile;
   cmfile.open("atoms-coords.meta");

   // check for open file
   if(!cmfile.is_open()){
      std::cerr << "Error! Coordinate metadata file atoms-coords.meta cannot be opened. Exiting" << std::endl;
      exit(1);
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
      ss >> vdc::num_atoms; // interpret as uint64_t
   }

   // Tell the user about the total number of atoms if verbosity is required
   if(vdc::verbose) std::cout << "   Number of atoms: " << num_atoms << std::endl;

   // get comment line
   getline(cmfile, line);

   //get number of materials
   int num_materials = 0;
   getline(cmfile, line);
   line.erase (line.begin(), line.begin()+20);
   num_materials=atoi(line.c_str());

   // Tell the user about the total number of materials if verbosity is required
   if(vdc::verbose) std::cout << "   Number of materials: " << num_materials << std::endl;

   // now read in material data
   for(int m = 0; m < num_materials; m++){
      vdc::material_t tmp; // temporary variable for reading data
      getline(cmfile, line);
      // convert to stringstream for reading
      std::istringstream ss(line);
      ss >> tmp.id >> tmp.moment >> tmp.element >> tmp.name;
      // save data to materials array
      vdc::materials.push_back(tmp);
      if(vdc::verbose) std::cout << "      Material " << tmp.id << "\t" << tmp.moment << "\t" << tmp.element << "\t" << tmp.name << std::endl;
   }

   // get comment line
   getline(cmfile, line);

   // get number of subsidiary files
   unsigned int num_coord_files;
   getline(cmfile, line);
   line.erase (line.begin(), line.begin()+22);
   num_coord_files=atoi(line.c_str());

   if(vdc::verbose) std::cout << "   Number of files: " << num_coord_files << std::endl;

   for(int file = 0; file < num_coord_files; file++){
      getline(cmfile, line);
      line.erase(remove(line.begin(), line.end(), '\t'), line.end());
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      line.erase(remove(line.begin(), line.end(), '\r'), line.end());
      vdc::coord_filenames.push_back(line);
      if(vdc::verbose) std::cout << "      " << line << std::endl;
   }

   return;

}

//------------------------------------------------------------------------------
// Function to read in coordinate data from subsidiary files
//------------------------------------------------------------------------------
void read_coord_data(){

   if(vdc::verbose) std::cout << "Reading coordinate data... " << std::flush;

   // resize arrays
   vdc::coordinates.resize(3*vdc::num_atoms);
   vdc::category.resize(vdc::num_atoms);
   vdc::type.resize(vdc::num_atoms);

   // index counter
   uint64_t atom_id = 0;

   // loop over all files
   for(unsigned int f = 0; f < vdc::coord_filenames.size(); f++){

      switch (vdc::format){

         case vdc::binary:{
            uint64_t num_atoms_in_file = 0;
            // open file in binary mode
            std::ifstream ifile;
            ifile.open(coord_filenames[f].c_str(), std::ios::binary); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << "Error! Coordinate data file \"" << coord_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
               std::cerr << "Error code: " << std::strerror(errno) << std::endl;
               exit(1);
            }
            // read number of atoms
            ifile.read( (char*)&num_atoms_in_file,sizeof(uint64_t) );
            // read type array
            ifile.read((char*)&vdc::type[atom_id], sizeof(int)*num_atoms_in_file);
            // read category array
            ifile.read((char*)&vdc::category[atom_id], sizeof(int)*num_atoms_in_file);
            ifile.read((char*)&vdc::coordinates[3*atom_id], sizeof(double)*num_atoms_in_file*3);
            // increment counter
            atom_id += num_atoms_in_file;
            ifile.close();

            break;
         }

         case vdc::text:{
            // open file
            std::ifstream ifile;
            ifile.open(coord_filenames[f].c_str()); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << "Error! Coordinate data file \"" << coord_filenames[f] << "\" cannot be opened. Exiting" << std::endl;
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
               vdc::type[atom_id] = type_id;
               vdc::category[atom_id] = category_id;
               vdc::coordinates[3*atom_id + 0] = x;
               vdc::coordinates[3*atom_id + 1] = y;
               vdc::coordinates[3*atom_id + 2] = z;
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

   //---------------------------------------------------------------
   // calculate system extent and centre
   //---------------------------------------------------------------
   double min[3] = {1e20, 1e20, 1e20};
   double max[3] = {0.0, 0.0, 0.0};
   double ave[3] = {0.0, 0.0, 0.0};

   for(unsigned int atom = 0; atom < vdc::num_atoms; atom++){

      // temporary variables
      double x = vdc::coordinates[3*atom + 0];
      double y = vdc::coordinates[3*atom + 1];
      double z = vdc::coordinates[3*atom + 2];

      // add coordinates to running total
      ave[0] += x;
      ave[1] += y;
      ave[2] += z;

      // calculate min and max
      if(x > max[0]) max[0] = x;
      if(y > max[1]) max[1] = y;
      if(z > max[2]) max[2] = z;

      if(x < min[0]) min[0] = x;
      if(y < min[1]) min[1] = y;
      if(z < min[2]) min[2] = z;

   }

   // save system dimensions
   vdc::system_size[0] = max[0] - min[0];
   vdc::system_size[1] = max[1] - min[1];
   vdc::system_size[2] = max[2] - min[2];

   // save system centre
   vdc::system_centre[0] = ave[0]/double(vdc::num_atoms);
   vdc::system_centre[1] = ave[1]/double(vdc::num_atoms);
   vdc::system_centre[2] = ave[2]/double(vdc::num_atoms);

   return;

}

}
