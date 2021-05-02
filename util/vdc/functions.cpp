//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Daniel Meilak 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include<cmath>
#include<algorithm>

// vdc headers
#include"vdc.hpp"

namespace vdc {

// forward function declaration
void arg_count(const input_t &input, size_t args_required, std::string requirement);

//------------------------------------------------------------------------------
// User specified initial frame to render
//------------------------------------------------------------------------------
void set_frame_start(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"frame-start\"\n\nExpects 1 argument: unsigned int\n\n"
                << "Sets initial frame for Povray rendering.\n\n"
                << "Default: frame-start = " << vdc::start_file_id << std::endl;
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "eq");

   int frame = std::stoi(input.value[0]);
   if (frame >= 0){ vdc::vdc_start_file_id = frame; }
   else {
      std::cerr << "Error - frame index cannot be negative." << std::endl;
      std::exit(EXIT_FAILURE);
   }
}


//------------------------------------------------------------------------------
// User specified final frame to render
//------------------------------------------------------------------------------
void set_frame_final(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"frame-final\"\n\nExpects 1 argument: unsigned int\n\n"
                << "Sets final frame for Povray rendering.\n"
                << "Default: frame-final = " << vdc::final_file_id << std::endl;
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "eq");

   int frame = std::stoi(input.value[0]);
   if (frame >= 0){ vdc::vdc_final_file_id = frame; }
   else {
      std::cerr << "Error - frame index cannot be negative." << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

//------------------------------------------------------------------------------
// User specified materials to remove
//------------------------------------------------------------------------------
void set_remove_materials(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"remove-materials\"\n\nExpects 1 or more arguments: unsigned int\n\n"
                << "Remove material IDs from visualisation. IDs start from 1.\n"
                << "Default: remove-materials = [empty]\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "ge");

   for (const std::string &mat : input.value){
      vdc::remove_materials.push_back(std::stoi(mat));
   }
}

//------------------------------------------------------------------------------
// User specified antiferromagnetic material spins to flip in povray visualisation
//------------------------------------------------------------------------------
void set_afm(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"afm\"\n\nExpects 1 or more arguments: unsigned int\n\n"
                << "Chosen material IDs are flipped in Povray visualisation. IDs start from 1.\n"
                << "Default: afm = [empty]\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "ge");

   for (const std::string &mat : input.value){
      vdc::afm_materials.push_back(std::stoi(mat));
   }
}

//------------------------------------------------------------------------
// Slice parameters
//------------------------------------------------------------------------
void set_slice(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice\"\n\nExpects 6 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined cube [xmin,xmax,ymin,ymax,zmin,zmax].\n"
                << "Default: slice = 0.0,1.0,0.0,1.0,0.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 6, "eq");

   // conver to double and store
   for (int i=0; i<6; i++){

      double val = std::stod(input.value[i]);
      if (val >= -0.000001 && val <= 1.000001){ vdc::slice_parameters[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (0,1) in '" << input.key << "' on line "
             << input.line_number << " of input file '" << vdc::input_file << std::endl;
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice;
}

//------------------------------------------------------------------------
// Slice void parameters
//------------------------------------------------------------------------
void set_slice_void(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-void\"\n\nExpects 6 arguments: real, in range (0,1)\n\n"
                << "Atoms inside defined cube [xmin,xmax,ymin,ymax,zmin,zmax] are removed.\n"
                << "Default: slice = 0.0,1.0,0.0,1.0,0.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 6, "eq");

   // conver to double and store
   for (int i=0; i<6; i++){

      double val = std::stod(input.value[i]);
      if (val >= -0.000001 && val <= 1.000001){ vdc::slice_parameters[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (0,1) in '" << input.key << "' on line "
             << input.line_number << " of input file '" << vdc::input_file << std::endl;
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_void;
}

//------------------------------------------------------------------------
// Slice sphere parameters
//------------------------------------------------------------------------
void set_slice_sphere(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-sphere\"\n\nExpects 3 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined ellipsoid [xfrac,yfrac,zfrac].\n"
                << "Default: slice = 1.0,1.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 3, "eq");

   // conver to double and store
   for (int i=0; i<3; i++){

      double val = std::stod(input.value[i]);
      if (val >= -0.000001 && val <= 1.000001){ vdc::slice_parameters[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (0,1) in '" << input.key << "' on line "
             << input.line_number << " of input file '" << vdc::input_file << std::endl;
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_sphere;
}

//------------------------------------------------------------------------
// Slice sphere parameters
//------------------------------------------------------------------------
void set_slice_cylinder(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-cylinder\"\n\nExpects 4 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined cylinder [xfrac,yfrac,zmin,zmax] along z-axis.\n"
                << "Default: slice = 1.0,1.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 4, "eq");

   // conver to double and store
   for (int i=0; i<4; i++){

      double val = std::stod(input.value[i]);
      if (val >= -0.000001 && val <= 1.000001){ vdc::slice_parameters[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (0,1) in '" << input.key << "' on line "
             << input.line_number << " of input file '" << vdc::input_file << std::endl;
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_cylinder;
}

void set_vector_z(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"vector-z\"\n\nExpects 3 arguments: real (brackets and comma optional)\n\n"
                << "Redefine z-axis orientation in cartesian coords.\n"
                << "Default: vector-z = {0.0,0.0,1.0}\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,3,"eq");

   for (int i=0; i<3; i++){
      vdc::vector_z[i] = std::stod(input.value[i]);
   }

   double length = std::sqrt(vector_z[0]*vector_z[0] + vector_z[1]*vector_z[1] + vector_z[2]*vector_z[2]);
   for (double &i : vector_z){
      i /= length;
   }

   vdc::z_vector = true;
}

void set_vector_x(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"vector-x\"\n\nExpects 3 arguments: real (brackets and comma optional)\n\n"
                << "Redefine x-axis orientation in cartesian coords. vector-z must also be defined.\n"
                << "Default: vector-x = {1.0,0.0,0.0}\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,3,"eq");

   for (int i=0; i<3; i++){
      vdc::vector_x[i] = std::stod(input.value[i]);
   }

   double length = std::sqrt(vector_x[0]*vector_x[0] + vector_x[1]*vector_x[1] + vector_x[2]*vector_x[2]);
   for (double &i : vector_x){
      i /= length;
   }

   vdc::x_vector = true;
}

void set_colourmap(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"colourmap\"\n\nExpects 1 argument: string\n\n"
                << "Choose preset colourmap. Availible maps: ";
                for (const std::string &map : vdc::colourmaps){
                   std::cout << map << " ";
                }
      std::cout << "\nDefault: colourmap = CBWR\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   if (std::find(vdc::colourmaps.begin(),vdc::colourmaps.end(),input.value[0]) != vdc::colourmaps.end() ){
      vdc::colour_keyword = input.value[0];
   }
   else {
      std::cerr << "Error - Colourmap keyword does not match." << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

void set_3D(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"3D\"\n\nExpects 1 argument: true or false\n\n"
                << "Enables brighness variation in Povray acorrding to spin component along vector-x.\n"
                << "Default: 3D = false\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   const std::string &value = input.value[0];

   if      (value == "true" ){ vdc::x_axis_colour = true; }
   else if (value == "false"){ vdc::x_axis_colour = false;}
   else {
      std::cerr << "Error - Expected true/false instead of " << value << " in " << input.key
                << "' on line " << input.line_number << " of input file '" << vdc::input_file << "'.\n";
      std::exit(EXIT_FAILURE);
   }
}

void set_custom_colourmap(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"custom-colourmap\"\n\nExpects 1 argument: filename\n\n"
                << "Povray uses User defined colourmap. Filename must contain 256 colours in RGB format:\n"
                << "3 columns of space separated reals in range (0,1).\n\n"
                << "Default: [empty]\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   // set colour_keyword to "custom"
   vdc::colour_keyword = "custom";

   // set custom map file name
   vdc::custom_colourmap_file = input.value[0];
}

// check number of args provided is correct
void arg_count(const input_t &input, size_t args_required, std::string requirement){

   const size_t &num_args = input.value.size();
   std::string many_or_few;

   if (requirement == "eq"){
      if (num_args == args_required){ return; }
      else {
         std::cerr << "Error - expected " << args_required << " arguments in " << input.key
                   << "' on line " << input.line_number << " of input file '" << vdc::input_file
                   << ".\nInstead got " << num_args << ".\n"; 
      }
   }
   else if (requirement == "ge"){
      if (num_args >= args_required){ return; }
      else {
         std::cerr << "Error - too few arguments passed to '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << std::endl;
      }
   }
   else if (requirement == "le"){
      if (num_args <= args_required){ return; }
      else {
         std::cerr << "Error - too many arguments passed to '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << std::endl;
      }
   }
   else {
      std::cerr << "Bad requirement: " << requirement << std::endl;
      std::exit(EXIT_FAILURE);
   }
   
   std::exit(EXIT_FAILURE);
}

} // end of namespace vdc