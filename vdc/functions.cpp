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
void error_message(const input_t &input, std::string message);

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

   // cannot be negative
   if (input.value[0][0] == '-'){ error_message(input,"frame index cannot be negative"); }

   // catch invalid conversion
   try { vdc::vdc_start_file_id = std::stoi(input.value[0]); }
   catch(...){ error_message(input,"invalid argument"); }

   // check frame range
   if (vdc::vdc_start_file_id > vdc::vdc_final_file_id){
      error_message(input,"start frame must be lower than start frame,");
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

   // cannot be negative
   if (input.value[0][0] == '-'){ error_message(input,"frame index cannot be negative"); }

   // catch invalid conversion
   try{ vdc::vdc_final_file_id = std::stoi(input.value[0]); }
   catch(...){ error_message(input,"invalid argument"); }

   // check frame range
   if (vdc::vdc_final_file_id < vdc::vdc_start_file_id){
      error_message(input,"final frame must be higher than start frame,");
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

   // add to remove_materials list
   for (const std::string &mat : input.value){

      // catch invalid conversion
      try{ vdc::remove_materials.push_back(std::stoi(mat)); }
      catch(...){ error_message(input,"invalid argument"); }
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

   // add to afm_materials list
   for (const std::string &mat : input.value){

      // catch invalid conversion
      try{ vdc::afm_materials.push_back(std::stoi(mat)); }
      catch(...){ error_message(input,"invalid argument"); }
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

   slice_t slice;
   slice.type = vdc::box;
   slice.param.resize(6);

   // convert to double and store
   for (int i=0; i<6; i++){

      // catch invalid conversion
      try { slice.param[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (slice.param[i] < -0.000001 || slice.param[i] > 1.000001){
         error_message(input,"fractional coords must be in range (0,1)");
      }
   }

   vdc::slices.push_back(slice);
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

   slice_t slice;
   slice.type = vdc::box_void;
   slice.param.resize(6);

   // conver to double and store
   for (int i=0; i<6; i++){

      // catch invalid conversion
      try { slice.param[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (slice.param[i] < -0.000001 || slice.param[i] > 1.000001){
         error_message(input,"fractional coords must be in range (0,1)");
      }
   }

   vdc::slices.push_back(slice);
}

//--------------------------------------------------------------------------------
// Slice sphere parameters cut system using x,y,z-axis coords, keeping inner atoms
//--------------------------------------------------------------------------------
void set_slice_sphere(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"slice-sphere\"\tExpects 3 arguments: positive real\n\n"
                << "Removes atoms inside of defined ellipsoid [xfrac,yfrac,zfrac].\n"
                << "Setting arguments to 1 excludes corners of system.\n\n"
                << "Default: [not set]\n"
                << "Example usage: slice-sphere = 0.5,1.0,1.0\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input, 3, "eq");

   slice_t slice;
   slice.type = vdc::sphere;
   slice.param.resize(3);

   // conver to double and store
   for (int i=0; i<3; i++){

      // catch invalid conversion
      try { slice.param[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (slice.param[i] < -0.000001){
         error_message(input,"fractional coords must be greater than 0");
      }
   }

   vdc::slices.push_back(slice);
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

   slice_t slice;
   slice.type = vdc::cylinder;
   slice.param.resize(4);

   // conver to double and store
   for (int i=0; i<4; i++){

      // catch invalid conversion
      try { slice.param[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (slice.param[i] < -0.000001 || slice.param[i] > 1.000001){
         error_message(input,"fractional coords must be in range (0,1)");
      }
   }

   vdc::slices.push_back(slice);
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
      try { vdc::vector_z[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }
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
      try { vdc::vector_x[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }
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
   else { error_message(input,"colourmap keyword does not match,"); }
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
   else { error_message(input,"expected true/false instead of '"+value+"'"); }
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

      // catch invalid conversion
      try{ vdc::camera_pos[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check argument range
      if (vdc::camera_pos[i] < -1.000001 || vdc::camera_pos[i] > 1.000001){
         error_message(input,"fractional coords must be in range (-1,1)");
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

      // catch invalid conversion
      try{ vdc::camera_look_at[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (vdc::camera_look_at[i] < -1.000001 || vdc::camera_look_at[i] > 1.000001){
         error_message(input,"fractional coords must be in range (-1,1)");
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
   try { vdc::camera_zoom = std::stod(input.value[0]); }
   catch(...){ error_message(input,"invalid argument"); }

   if (vdc::camera_zoom <= 0.0){
      error_message(input,"camera zoom must be greater than 0.0,");
   }
}

//----------------------------------------------------------------------------------
// Set povray stick cutoff
//----------------------------------------------------------------------------------
void set_sticks_cutoff(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"sticks-cutoff\"\tExpects 1 argument: positive real\n\n"
                << "Povray sticks cutoff. Sticks connect nearby atoms to better\n"
                << "visualise the atomic structure in povray. The cutoff gives\n"
                << "the maximum distance in Angstroms between atoms connected\n"
                << "sticks. Note that this option is very inefficient for large\n"
                << "system sizes."
                << "\nExample usage: sticks-cutoff = 3.54\t\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,1,"eq");

   // convert to double and store
   vdc::sticks_cutoff = std::stod(input.value[0]);

   // check range
   if (vdc::sticks_cutoff <= 0.0){
      error_message(input,"sticks-cutoff must be greater than 0.0");
   }

   return;

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

      // catch invalid conversion
      try { vdc::atom_sizes.push_back(std::stod(value)); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (vdc::atom_sizes.back() <= 0.0){
         error_message(input,"atom size must be greater than 0.0");
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

      // catch invalid conversion
      try { vdc::arrow_sizes.push_back(std::stod(value)); }
      catch(...){ error_message(input,"invalid argument"); }

      // check range
      if (vdc::arrow_sizes.back() <= 0.0){
         error_message(input,"atom size must be greater than 0.0");
      }
   }
}

//----------------------------------------------------------------------------------
// Set cell size components for cells output
//----------------------------------------------------------------------------------
void set_cell_size(const input_t &input){

   // print help message if argument is "-h"
   if (input.value[0] == "-h"){
      std::cout << "\"cell-size\"\tExpects 3 arguments: real (brackets and comma optional)\n\n"
                << "Set cell sizes in x,y,z in Angstroms\n\n"
                << "Default: cell-size = {10.0,10.0,10.0}\n";
      std::exit(EXIT_SUCCESS);
   }

   // check args
   arg_count(input,3,"eq");

   for (int i=0; i<3; i++){
      try { vdc::cell_size[i] = std::stod(input.value[i]); }
      catch(...){ error_message(input,"invalid argument"); }
   }

}

// bookkeeping functino to check number of args provided is correct
void arg_count(const input_t &input, size_t args_required, std::string requirement){

   const size_t &num_args = input.value.size();
   std::string many_or_few;

   // location message for input file or command line parameters
   std::string location_message;
   if (input.line_number == -1){ location_message = "' passed in the command line\n"; }
   else { location_message = "' on line " + std::to_string(input.line_number) + " of input file '" + vdc::input_file + "'\n"; }

   if (requirement == "eq"){
      if (num_args == args_required){ return; }
      else {
         std::cerr << "Error - expected " << args_required << " arguments in '" << input.key
                   << location_message << "Instead got " << num_args << "\n";
      }
   }
   else if (requirement == "ge"){
      if (num_args >= args_required){ return; }
      else { std::cerr << "Error - too few arguments passed to '" << input.key << location_message; }
   }
   else if (requirement == "le"){
      if (num_args <= args_required){ return; }
      else { std::cerr << "Error - too many arguments passed to '" << input.key << location_message; }
   }
   else {
      std::cerr << "Bad requirement: " << requirement << std::endl;
      std::exit(EXIT_FAILURE);
   }

   std::exit(EXIT_FAILURE);
}

// output error message and exit
void error_message(const input_t &input, std::string message){

   // location message for input file or command line parameters
   std::string location_message;
   if (input.line_number == -1){ location_message = "' passed in the command line\n"; }
   else { location_message = "' on line " + std::to_string(input.line_number) + " of input file '" + vdc::input_file + "'\n"; }

   std::cerr << "Error - " << message << " in '" << input.key << location_message;

   std::exit(EXIT_FAILURE);
}

} // end of namespace vdc
