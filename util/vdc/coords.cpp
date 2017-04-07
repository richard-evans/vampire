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
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
void read_coord_metadata(std::vector <std::string>& filenames);

//------------------------------------------------------------------------------
// Wrapper function to read coordinate metafile to initialise data structures
// and process coordinate data
//------------------------------------------------------------------------------
void process_coordinates(){

   // array to store subsidiary data file names
   std::vector <std::string> filenames(0);

   // read meta data
   vdc::read_coord_metadata(filenames);

   // load coordinates into buffers

   vdc::coordinates.resize(3*vdc::num_atoms);
   vdc::category.resize(vdc::num_atoms);
   vdc::type.resize(vdc::num_atoms);

   // index counter
   uint64_t atom_id = 0;

   // loop over all files
   for(unsigned int f = 0; f < filenames.size(); f++){

      switch (vdc::format){

         case vdc::binary:
            break;

         case vdc::text:{
            // open file
            std::ifstream ifile;
            ifile.open(filenames[f].c_str()); // check for errors
            // check for open file
            if(!ifile.is_open()){
               std::cerr << "Error! Coordinate data file \"" << filenames[f] << "\" cannot be opened. Exiting" << std::endl;
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

   // output xyz file
   std::ofstream ofile;
   ofile.open("crystal.xyz");
   ofile << vdc::num_atoms << "\n\n";

   for(uint64_t atom = 0; atom < vdc::num_atoms; atom++){
      int type_id = vdc::type[atom];
      ofile << materials[type_id].element << "\t" << vdc::coordinates[3*atom + 0] << "\t" << vdc::coordinates[3*atom + 1] << "\t" << vdc::coordinates[3*atom + 2] << "\n";
   }

   ofile << std::flush;
   ofile.close();

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
void read_coord_metadata(std::vector <std::string>& filenames){

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
   if(data_format_str == test) vdc::format = vdc::text;
   test = "binary";
   if(data_format_str == test) vdc::format = vdc::binary;
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

   // Tell the user about the total number of atoms if verbosity is required
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
      filenames.push_back(line);
      if(vdc::verbose) std::cout << "      " << line << std::endl;
   }

   return;

}

}
