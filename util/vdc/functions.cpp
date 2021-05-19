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
      std::cout << "\"frame-start\"\tExpects 1 argument: unsigned int\n\n"
                << "Sets initial frame for Povray rendering.\n\n"
                << "Default: frame-start = " << vdc::start_file_id << std::endl;
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "eq");

   int frame = std::stoi(input.value[0]);
   if (frame >= 0){ vdc::vdc_start_file_id = frame; }
   else {
      std::cerr << "Error - frame index cannot be negative in '" << input.key << "' on line "
                << input.line_number << "of input file '" << vdc::input_file << "'\n";
      std::exit(EXIT_FAILURE);
   }
}


//------------------------------------------------------------------------------
// User specified final frame to render
//------------------------------------------------------------------------------
void set_frame_final(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"frame-final\"\tExpects 1 argument: unsigned int\n\n"
                << "Sets final frame for Povray rendering.\n\n"
                << "Default: frame-final = " << vdc::final_file_id << std::endl;
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "eq");

   int frame = std::stoi(input.value[0]);
   if (frame >= 0){ vdc::vdc_final_file_id = frame; }
   else {
      std::cerr << "Error - frame index cannot be negative in '" << input.key << "' on line "
                << input.line_number << "of input file '" << vdc::input_file << "'\n";
      std::exit(EXIT_FAILURE);
   }
}

//------------------------------------------------------------------------------
// User specified materials to remove
//------------------------------------------------------------------------------
void set_remove_materials(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"remove-materials\"\tExpects 1 or more arguments: unsigned int\n\n"
                << "Remove material IDs from visualisation. IDs start from 1.\n\n"
                << "Default: [not set]\n"
                << "Example usage: remove-materials = 1,2\n";
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
      std::cout << "\"afm\"\tExpects 1 or more arguments: unsigned int\n\n"
                << "Chosen material IDs are flipped in Povray visualisation. IDs start from 1.\n\n"
                << "Default: [not set]\n"
                << "Example usage: afm = 1,2\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 1, "ge");

   for (const std::string &mat : input.value){
      vdc::afm_materials.push_back(std::stoi(mat));
   }
}

//--------------------------------------------------------------------------------
// Slice parameters cut system along x,y,z-axis directions, keeping inner atoms
//--------------------------------------------------------------------------------
void set_slice(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice\"\t\tExpects 6 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined cube [xmin,xmax,ymin,ymax,zmin,zmax].\n\n"
                << "Default: slice = 0.0,1.0,0.0,1.0,0.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 6, "eq");

   // convert to double and store
   for (int i=0; i<6; i++){

      double val = std::stod(input.value[i]);
      if (val >= -0.000001 && val <= 1.000001){ vdc::slice_parameters[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (0,1) in '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice;
}

//--------------------------------------------------------------------------------------
// Slice void parameters cut system along x,y,z-axis directions, keeping outer atoms
//--------------------------------------------------------------------------------------
void set_slice_void(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-void\"\tExpects 6 arguments: real, in range (0,1)\n\n"
                << "Atoms inside defined cube [xmin,xmax,ymin,ymax,zmin,zmax] are removed.\n\n"
                << "Default: [not set]\n"
                << "Example usage: slice-void = 0.1,0.9,0.1,0.9,0.1,0.9\n";
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
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_void;
}

//--------------------------------------------------------------------------------
// Slice sphere parameters cut system using x,y,z-axis coords, keeping inner atoms
//--------------------------------------------------------------------------------
void set_slice_sphere(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-sphere\"\tExpects 3 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined ellipsoid [xfrac,yfrac,zfrac].\n\n"
                << "Default: [not set]\n"
                << "Example usage: slice-sphere = 0.5,1.0,1.0\n";
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
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_sphere;
}

//----------------------------------------------------------------------------------
// Slice cylinder parameters cut system using x,y,z-axis coords, keeping inner atoms
//----------------------------------------------------------------------------------
void set_slice_cylinder(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-cylinder\"\tExpects 4 arguments: real, in range (0,1)\n\n"
                << "Removes atoms outside of defined cylinder [xfrac,yfrac,zmin,zmax] along z-axis.\n\n"
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
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }

   vdc::slice_type = vdc::slice_cylinder;
}

//----------------------------------------------------------------------------------
// Set z-vector direction which defines principal colour axis
//----------------------------------------------------------------------------------
void set_vector_z(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"vector-z\"\tExpects 3 arguments: real (brackets and comma optional)\n\n"
                << "Redefine z-axis orientation in cartesian coords.\n\n"
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

//----------------------------------------------------------------------------------
// Set x-vector direction which defines secondary colour axis
//----------------------------------------------------------------------------------
void set_vector_x(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"vector-x\"\tExpects 3 arguments: real (brackets and comma optional)\n\n"
                << "Redefine x-axis orientation in cartesian coords. vector-z must also be defined.\n\n"
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

//----------------------------------------------------------------------------------
// Choose predefined colourmap
//----------------------------------------------------------------------------------
void set_colourmap(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"colourmap\"\tExpects 1 argument: string\n\n"
                << "Choose preset colourmap. Availible maps: ";
                for (const std::string &map : vdc::colourmaps){
                   std::cout << map << " ";
                }
      std::cout << "\n\nDefault: colourmap = CBWR\n";
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

//----------------------------------------------------------------------------------
// Use 3D colouring scheme
//----------------------------------------------------------------------------------
void set_3D(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"3D\"\tExpects 1 argument: true or false\n\n"
                << "Enables brighness variation in Povray acorrding to spin component along vector-x.\n\n"
                << "Default: [not set]\n"
                << "Example usage: 3D = true\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   const std::string &value = input.value[0];

   if      (value == "true" ){ vdc::x_axis_colour = true; }
   else if (value == "false"){ vdc::x_axis_colour = false;}
   else {
      std::cerr << "Error - Expected true/false instead of " << value << " in " << input.key
                << "' on line " << input.line_number << " of input file '" << vdc::input_file << "'\n";
      std::exit(EXIT_FAILURE);
   }
}

//----------------------------------------------------------------------------------
// Set user-provided colourmap
//----------------------------------------------------------------------------------
void set_custom_colourmap(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"custom-colourmap\"\tExpects 1 argument: filename\n\n"
                << "Povray uses User defined colourmap. Filename must contain 256 colours in RGB format:\n"
                << "3 columns of space separated reals in range (0,1).\n\n"
                << "Default: [not set]\n"
                << "Example usage: custom_colourmap = colourmapfile.txt \n";       
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   // set colour_keyword to "custom"
   vdc::colour_keyword = "custom";

   // set custom map file name
   vdc::custom_colourmap_file = input.value[0];
}

//----------------------------------------------------------------------------------
// Set povray camera location
//----------------------------------------------------------------------------------
void set_camera_position(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"camera-position\"\tExpects 3 arguments: real, in range (-1,1)\n\n"
                << "Povray camera position, set using fractional coords.\n"
                << "Camera distance from look at point is calculated automatically however\n"
                << "it can be changed by using camera-zoom.\n\n"
                << "Default: [not set]\n"
                << "Example usage: camera-position = (1,1,1)\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,3,"eq");

   // convert to double and store
   for (int i=0; i<3; i++){

      double val = std::stod(input.value[i]);
      if (val >= -1.000001 && val <= 1.000001){ vdc::camera_pos[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (-1,1) in '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }

   // set default camera pos flag to false
   vdc::default_camera_pos = false;
}

//----------------------------------------------------------------------------------
// Set povray camera look at position
//----------------------------------------------------------------------------------
void set_camera_look_at(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"camera-look-at\"\tExpects 3 arguments: real, in range (-1,1)\n\n"
                << "Povray camera look at position, set using fractional coords.\n"
                << "The position is a location in the bounding box of the system, with centre (0,0,0).\n\n"
                << "Default: camera-look-at = (0.0,0.0,0.0)\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,3,"eq");

   // convert to double and store
   for (int i=0; i<3; i++){

      double val = std::stod(input.value[i]);
      if (val >= -1.000001 && val <= 1.000001){ vdc::camera_look_at[i] = val; }
      else {
         std::cerr << "Error - fractional coords must be in range (-1,1) in '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }
}

//----------------------------------------------------------------------------------
// Set povray camera distance multiplier
//----------------------------------------------------------------------------------
void set_camera_zoom(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"camera-zoom\"\tExpects 1 argument: positive real\n\n"
                << "Povray camera zoom multiplier.\n"
                << "The default distance from the camera is automatically calculated according to the size of the system.\n"
                << "This can be increased or reduced using camera-zoom to multiply the default distance. Values lower\n"
                << "than 1.0 reduce the distance, while values above 1.0 increase it.\n\n"
                << "Default: camera-zoom = 1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   // set value, must be greater than 0
   vdc::camera_zoom = std::stod(input.value[0]);
   if (vdc::camera_zoom <= 0.0){
      std::cerr << "Error - camera zoom must be greater than 0.0, in '" << input.key << "' on line "
                << input.line_number << "of input file '" << vdc::input_file << "'\n";
      std::exit(EXIT_FAILURE);
   }
}

//----------------------------------------------------------------------------------
// Set povray background colour
//----------------------------------------------------------------------------------
void set_background_colour(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"background-colour\"\tExpects 1 argument: positive real\n\n"
                << "Povray background colour.\n"
                << "Povray includes various predefined colours such as \"White\",\"Black\",\"Gray\".\n"
                << "A list of these names and their rgb values can be found at: https://github.com/POV-Ray/povray/blob/master/distribution/include/colors.inc\n"
                << "Note: misspelled colour names will not be detected by vdc but will cause errors in Povray.\n\n"
                << "Default: background-colour = " << vdc::background_colour << "\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   // set value
   vdc::background_colour = input.value[0];

   // preset colour in Povray are capitalised
   vdc::background_colour[0] = std::toupper(vdc::background_colour[0]);
}

//----------------------------------------------------------------------------------
// Set povray atom size
//----------------------------------------------------------------------------------
void set_atom_sizes(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"atom-sizes\"\tExpects 1 or more arguments: positive real\n\n"
                << "Povray atom sizes. Atoms are represented by spheres with a defined radius.\n"
                << "Individual materials can have different atom sizes.\n\n"
                << "Default: all atoms are set to " << vdc::atom_sizes[0] 
                << "\nExample usage: atom-sizes = {1.2,2.0}\t// set mat1: 1.2 and mat2: 2.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"ge");

   // convert to double and store
   vdc::atom_sizes.clear();
   for (const std::string &value : input.value){

      vdc::atom_sizes.push_back(std::stod(value));
      if (vdc::atom_sizes.back() <= 0.0){
         std::cerr << "Error - atom size must be greater than 0.0 in '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }
}

//----------------------------------------------------------------------------------
// Set povray arrow size
//----------------------------------------------------------------------------------
void set_arrow_sizes(const input_t &input){

   // print help message if argument is "-h"
      if (input.value[0] == "-h"){
      std::cout << "\"arrow-sizes\"\tExpects 1 or more arguments: positive real\n\n"
                << "Povray arrow sizes.\n"
                << "Individual materials can have different arrow sizes.\n\n"
                << "Default: all arrows are set to " << vdc::arrow_sizes[0] 
                << "\nExample usage: arrow-sizes = {1.2,2.0}\t// set mat1: 1.2 and mat2: 2.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"ge");

   // convert to double and store
   vdc::arrow_sizes.clear();
   for (const std::string &value : input.value){

      vdc::arrow_sizes.push_back(std::stod(value));
      if (vdc::arrow_sizes.back() <= 0.0){
         std::cerr << "Error - arrow size must be greater than 0.0 in '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
         std::exit(EXIT_FAILURE);  
      }
   }
}

// bookkeeping functino to check number of args provided is correct
void arg_count(const input_t &input, size_t args_required, std::string requirement){

   const size_t &num_args = input.value.size();
   std::string many_or_few;

   if (requirement == "eq"){
      if (num_args == args_required){ return; }
      else {
         std::cerr << "Error - expected " << args_required << " arguments in " << input.key
                   << "' on line " << input.line_number << " of input file '" << vdc::input_file
                   << "'\nInstead got " << num_args << ".\n"; 
      }
   }
   else if (requirement == "ge"){
      if (num_args >= args_required){ return; }
      else {
         std::cerr << "Error - too few arguments passed to '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
      }
   }
   else if (requirement == "le"){
      if (num_args <= args_required){ return; }
      else {
         std::cerr << "Error - too many arguments passed to '" << input.key << "' on line "
                   << input.line_number << " of input file '" << vdc::input_file << "'\n";
      }
   }
   else {
      std::cerr << "Bad requirement: " << requirement << std::endl;
      std::exit(EXIT_FAILURE);
   }
   
   std::exit(EXIT_FAILURE);
}

} // end of namespace vdc