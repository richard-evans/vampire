//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans, Daniel Meilak 2017-2019. All rights reserved.
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
void calculate_system_extent(std::vector<int>& magnetic_list, std::vector<int>& non_magnetic_list);
void slice_system();
bool box_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound);
bool sphere_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound);
bool cylinder_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound);


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
   vdc::read_nm_metadata();
   vdc::read_nm_data();

   // Calculate system dimensions
   vdc::calculate_system_extent(vdc::atoms_list,vdc::nm_atoms_list);

   // Create list of atoms in slice
   vdc::slice_system();

   // Calculate systenm dimensions after slicing
   vdc::calculate_system_extent(vdc::sliced_atoms_list,vdc::sliced_nm_atoms_list);

   // output xyz file;
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

   for(unsigned int file = 0; file < num_coord_files; file++){
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

   // create vector list of atom indices
   vdc::atoms_list.resize(vdc::num_atoms);
   for(unsigned int atom = 0; atom < vdc::num_atoms; atom++){
      vdc::atoms_list.push_back(atom);
   }

   return;
}

//------------------------------------------------------------------------------
// calculate system extent and centre
//------------------------------------------------------------------------------
void calculate_system_extent(std::vector<int>& magnetic_list, std::vector<int>& non_magnetic_list){

   double min[3] = {1e20, 1e20, 1e20};
   double max[3] = {0.0, 0.0, 0.0};
   //double ave[3] = {0.0, 0.0, 0.0};

   // loop through all magnetic atoms
   for(size_t i=0; i < magnetic_list.size(); i++){

      // get atom ID
      unsigned int atom = magnetic_list[i];

      // temporary variables
      double x = vdc::coordinates[3*atom + 0];
      double y = vdc::coordinates[3*atom + 1];
      double z = vdc::coordinates[3*atom + 2];

      // calculate min and max
      if(x > max[0]) max[0] = x;
      if(y > max[1]) max[1] = y;
      if(z > max[2]) max[2] = z;

      if(x < min[0]) min[0] = x;
      if(y < min[1]) min[1] = y;
      if(z < min[2]) min[2] = z;
   }

   // loop through all non-magnetic atoms
   for(size_t i=0; i < non_magnetic_list.size(); i++){

      // get atom ID
      unsigned int atom = non_magnetic_list[i];

      // temporary variables
      double x = vdc::nm_coordinates[3*atom + 0];
      double y = vdc::nm_coordinates[3*atom + 1];
      double z = vdc::nm_coordinates[3*atom + 2];

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
   vdc::system_centre[0] = (max[0] + min[0])/2.0;
   vdc::system_centre[1] = (max[1] + min[1])/2.0;
   vdc::system_centre[2] = (max[2] + min[2])/2.0;

   return;
}

//---------------------------------------------------------------
// Find list of atoms in user defined slice
//---------------------------------------------------------------
void slice_system(){
   
   // work out borders for slice param
   for (slice_t &slice : vdc::slices){

      switch (slice.type){
      // group box and box void as min and max ranges are the same
      case vdc::box :
      case vdc::box_void :
         slice.bound.resize(6); // xmin,xmax,ymin,ymax,zmin,zmax

         // min
         slice.bound[0] = (slice.param[0]*vdc::system_size[0])-(vdc::system_size[0]*0.5)+vdc::system_centre[0];
         slice.bound[2] = (slice.param[2]*vdc::system_size[1])-(vdc::system_size[1]*0.5)+vdc::system_centre[1];
         slice.bound[4] = (slice.param[4]*vdc::system_size[2])-(vdc::system_size[2]*0.5)+vdc::system_centre[2];

         // max
         slice.bound[1] = (slice.param[1]*vdc::system_size[0])-(vdc::system_size[0]*0.5)+vdc::system_centre[0];
         slice.bound[3] = (slice.param[3]*vdc::system_size[1])-(vdc::system_size[1]*0.5)+vdc::system_centre[1];
         slice.bound[5] = (slice.param[5]*vdc::system_size[2])-(vdc::system_size[2]*0.5)+vdc::system_centre[2];
         break;
      
      case vdc::sphere :
         slice.bound.resize(3); // a,b,c

         // work out radii of the ellipse
         slice.bound[0] = vdc::system_size[0]*slice.param[0]/2.0;
         slice.bound[1] = vdc::system_size[1]*slice.param[1]/2.0;
         slice.bound[2] = vdc::system_size[2]*slice.param[2]/2.0;
         break;
      
      case vdc::cylinder :
         slice.bound.resize(4); // a,b,zmin,zmax

         slice.bound[0] = vdc::system_size[0]*slice.param[0]/2.0;
         slice.bound[1] = vdc::system_size[1]*slice.param[1]/2.0;
         slice.bound[2] = (slice.param[2]*vdc::system_size[2])-(vdc::system_size[2]*0.5)+vdc::system_centre[2];
         slice.bound[3] = (slice.param[3]*vdc::system_size[2])-(vdc::system_size[2]*0.5)+vdc::system_centre[2];
         break;
      }
   }

   // find all valid atoms in slices
   for (unsigned int atom = 0; atom < vdc::num_atoms; atom++){

      // remove materials defined by user
      if ( std::find(vdc::remove_materials.begin(), vdc::remove_materials.end(), vdc::type[atom]+1) != vdc::remove_materials.end() ){
         continue;
      }

      // if no slices defined, include all atoms
      if ( vdc::slices.empty() ){
         vdc::sliced_atoms_list.push_back(atom);
      }
      // otherwise include all atoms in any slice
      else {

         // atom coords
         double x = vdc::coordinates[3*atom + 0];
         double y = vdc::coordinates[3*atom + 1];
         double z = vdc::coordinates[3*atom + 2];

         for (slice_t slice : slices){

            // flag to show atom is in any slice
            bool in_bounds = false;

            switch (slice.type){
            case vdc::box :
               in_bounds = box_slice(x,y,z,slice.bound);
               break;
            
            case vdc::box_void :
               in_bounds = !box_slice(x,y,z,slice.bound);
               break;
            
            case vdc::sphere :
               in_bounds = sphere_slice(x,y,z,slice.bound);
               break;

            case vdc::cylinder :
               in_bounds = cylinder_slice(x,y,z,slice.bound);
               break;

            default:
               std::cerr << "Error - unknown slice type, probably a mistake in command.cpp" << std::endl;
               std::exit(EXIT_FAILURE);
            }

            // if the atom is in any slice, add to final list and stop checking others
            if (in_bounds){ 
               vdc::sliced_atoms_list.push_back(atom);
               break;
            }
         }
      }
   }

   // find all valid non_magnetic atoms in slices
   for (unsigned int atom = 0; atom < vdc::num_nm_atoms; atom++){

      // remove materials defined by user
      if ( std::find(vdc::remove_materials.begin(), vdc::remove_materials.end(), vdc::nm_type[atom]+1) != vdc::remove_materials.end() ){
         continue;
      }

      // if no slices defined, include all atoms
      if ( vdc::slices.empty() ){
         vdc::sliced_nm_atoms_list.push_back(atom);
      }
      // otherwise include all atoms in any slice
      else {

         // atom coords
         double x = vdc::nm_coordinates[3*atom + 0];
         double y = vdc::nm_coordinates[3*atom + 1];
         double z = vdc::nm_coordinates[3*atom + 2];

         for (slice_t slice : slices){

            // flag to show atom is in any slice
            bool in_bounds = false;

            switch (slice.type){
            case vdc::box :
               in_bounds = box_slice(x,y,z,slice.bound);
               break;
            
            case vdc::box_void :
               in_bounds = !box_slice(x,y,z,slice.bound);
               break;
            
            case vdc::sphere :
               in_bounds = sphere_slice(x,y,z,slice.bound);
               break;

            case vdc::cylinder :
               in_bounds = cylinder_slice(x,y,z,slice.bound);
               break;

            default:
               std::cerr << "Error - unknown slice type, probably a mistake in command.cpp" << std::endl;
               std::exit(EXIT_FAILURE);
            }

            // if the atom is in any slice, add to final list and stop checking others
            if (in_bounds){ 
               vdc::sliced_nm_atoms_list.push_back(atom);
               break;
            }
         }
      }
   } 

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;
}

bool box_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound){

   return (x >= bound[0] && x <= bound[1]) && (y >= bound[2] && y <= bound[3]) && (z >= bound[4] && z <= bound[5]);   
}

bool sphere_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound){

   double temp = 0;

   temp += (x - vdc::system_centre[0])*(x - vdc::system_centre[0])/(bound[0]*bound[0]);
   temp += (y - vdc::system_centre[1])*(y - vdc::system_centre[1])/(bound[1]*bound[1]);
   temp += (z - vdc::system_centre[2])*(z - vdc::system_centre[2])/(bound[2]*bound[2]);

   return temp <= 1.0;
}

bool cylinder_slice(const double &x, const double &y, const double &z, const std::vector<double> &bound){

   double temp = 0;

   temp += (x - vdc::system_centre[0])*(x - vdc::system_centre[0])/(bound[0]*bound[0]);
   temp += (y - vdc::system_centre[1])*(y - vdc::system_centre[1])/(bound[1]*bound[1]);

   return temp<=1.0 && z>=bound[2] && z<=bound[3];
}

} // end of namespace vdc
