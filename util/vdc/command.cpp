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
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
int extract_vector( std::string arg_string, std::vector<double>& arg_vector );
int extract_materials( std::string arg_string, std::vector<int>& arg_vector);
int extract_slice_param( std::string arg_string, std::vector<double>& arg_vector, int number_of_param);
int check_arg( int& arg, int argc, char* argv[], std::string& temp_str, std::string error_output );
void init_vector_y();

//------------------------------------------------------------------------------
// Command line parsing function
//------------------------------------------------------------------------------
int command( int argc, char* argv[] ){
// Command line options for utility to be implemented:
//    --xyz - generate xyz file
//    --povray - generate povray files
//    --vector - generate raw xyz vector data
//    --vtk - generate vtk files
//    --cells - collate data into cells
//    --cell-size = x - define cell size
//    --povray-cells - generate renderable cell positions for povray
//    --verbose [= true, false] - set verbose output to screen
//    --colours = default, rwb [red-white-blue], oyb [orange-yellow-blue], jet, gs [grey-scale], cw [colour-wheel], <filename>
//    --objects = spins, cones, spheres, cubes
//    --slice = x,x,y,y,z,z
//    --multiscale = gradient, material, region

   // temporary string for storing command line argument
   std::string temp_str;

   for (int arg = 1; arg < argc; arg++){

      // read prefix
      std::string sw=argv[arg];

      //------------------------------------------------------------------------
      // Check for appropriate data outputs
      //------------------------------------------------------------------------
      // xyz coordinate file output
      if (sw == "--xyz"){
         vdc::xyz = true;
      }
      // xyz coordinate file output
      else if (sw == "--povray"){
         vdc::povray = true;
      }
      // xyz coordinate file output
      else if (sw == "--vtk"){
         vdc::vtk = true;
      }
      // plain text file output
      else if (sw == "--text"){
         vdc::txt = true;
      }
      //------------------------------------------------------------------------
      // Check for verbose output
      //------------------------------------------------------------------------
      else if (sw == "--verbose"){
         vdc::verbose = true;
      }
      //------------------------------------------------------------------------
      // Check for user specified initial frame to render
      //------------------------------------------------------------------------
      else if (sw == "--frame-start"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected index of intial frame to render." );

         if ( stoi(temp_str) >= 0 ){
            vdc::vdc_start_file_id = stoi(temp_str);
         }
         else{
            ////terminaltextcolor(RED);
            std::cerr << "Error - frame index cannot be negative."
                      << std::endl;
            return EXIT_FAILURE;
            ////terminaltextcolor(WHITE);
         }
      }
      //------------------------------------------------------------------------
      // Check for user specified final frame to render
      //------------------------------------------------------------------------
      else if (sw == "--frame-final"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected index of final frame to render." );

         if ( stoi(temp_str) >= 0 ){
            vdc::vdc_final_file_id = stoi(temp_str)+1; // +1 since we want to render the last frame as well
         }
         else{
            ////terminaltextcolor(RED);
            std::cerr << "Error - frame index cannot be negative."
                      << std::endl;
            return EXIT_FAILURE;
            ////terminaltextcolor(WHITE);
         }
      }
      //------------------------------------------------------------------------
      // Check for user specified materials to remove
      //------------------------------------------------------------------------
      else if (sw == "--remove-material" || sw == "--remove-materials"){

            // check number of args not exceeded
            check_arg(arg, argc, argv, temp_str, "Error - expected at least one variable." );

            // work through vector and extract values
            extract_materials(temp_str, vdc::remove_materials);
      }
      //------------------------------------------------------------------------
      // Check for user specified antiferromagnetic material spins to flip in
      // povray visualisation
      //------------------------------------------------------------------------
      else if (sw == "--afm"){

            // check number of args not exceeded
            check_arg(arg, argc, argv, temp_str, "Error - expected at least one variable." );

            // work through vector and extract values
            extract_materials(temp_str, vdc::afm_materials);
      }
      //------------------------------------------------------------------------
      // Check for slice parameters
      //------------------------------------------------------------------------
      else if (sw == "--slice"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 6 comma separated variables." );

         // work through vector and extract values
         extract_slice_param(temp_str, vdc::slice_parameters, 6);

         // set slice keyword
         vdc::slice_type = "slice";
      }
      else if (sw == "--slice-void"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 6 comma separated variables." );

         // work through vector and extract values
         extract_slice_param(temp_str, vdc::slice_parameters, 6);

         // set slice keyword
         vdc::slice_type = "slice-void";
      }
      else if (sw == "--slice-sphere"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

         // work through vector and extract values
         extract_slice_param(temp_str, vdc::slice_parameters, 3);

         // set slice keyword
         vdc::slice_type = "slice-sphere";
      }
      else if (sw == "--slice-cylinder"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 4 comma separated variables." );

         // work through vector and extract values
         extract_slice_param(temp_str, vdc::slice_parameters, 4);

         // set slice keyword
         vdc::slice_type = "slice-cylinder";
      }
      //------------------------------------------------------------------------
      // Check for colour mapping parameters
      //------------------------------------------------------------------------
      else if (sw == "--vector-z"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

         // work through vector and extract values
         extract_vector(temp_str, vdc::vector_z );

         // confirm initialisation of z-axis
         z_vector = true;
      }
      else if (sw == "--vector-x"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

         // work through vector and extract values
         extract_vector(temp_str, vdc::vector_x );

         // confirm initialisation of x-axis
         x_vector = true;
      }
      else if (sw == "--colourmap"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected colourmap keyword." );

         if ( temp_str == "C2" ){
            vdc::colour_keyword = temp_str;
         }
         else if (temp_str == "BWR" ){
            vdc::colour_keyword = temp_str;
         }
         else if (temp_str == "CBWR" ){
            vdc::colour_keyword = temp_str;
         }
         else if (temp_str == "Rainbow" ){
            vdc::colour_keyword = temp_str;
         }
         else {
            ////terminaltextcolor(RED);
            std::cerr << "Error - Colourmap keyword does not match."
                      << std::endl;
            return EXIT_FAILURE;
            ////terminaltextcolor(WHITE);
         }

      }
      else if ( sw == "--3D" ){

         vdc::x_axis_colour = true;

      }
      else if ( sw == "--custom-colourmap" ){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected custom colourmap name.");

         // set colour_keyword to "custom" 
         vdc::colour_keyword = "custom";

         // set custom map file name
         vdc::custom_colourmap_file = temp_str;

      }
      else {
         ////terminaltextcolor(RED);
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'"
                   << std::endl;
         ////terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
   }

   //---------------------------------------------------------------------------
   // Additional checks on input parameters
   //---------------------------------------------------------------------------

   // check that some kind of data output is requested
   if( !vdc::xyz && !vdc::povray && !vdc::vtk && !vdc::txt ){
      std::cerr << "Error! No output data formats requested. Available options are: " << std::endl;
      std::cerr << "\t\t --xyz    Data output in .xyz format for viewing in rasmol/jmol" << std::endl;
      std::cerr << "\t\t --povray Data output in PoVRAY format for rendering" << std::endl;
      std::cerr << "\t\t --vtk    Data output in VTK format for viewing in Paraview" << std::endl;
      std::cerr << "\t\t --text   Data output in plain text format for plotting in gnuplot/excel etc" << std::endl;
      return EXIT_FAILURE;
   }

   // check for valid axis initialisations
   if ( z_vector && !x_vector ){

      // check for a z-axis with vector_z[2] = 0
      // (Hence z_vector lies in xy-plane)
      if ( (-0.000001 < vdc::vector_z[2]) && (vdc::vector_z[2] < 0.000001) ){
         // x-axis will lie along {0,0,1}
         vdc::vector_x = {0.0, 0.0, 1.0};

         // find vector_y
         init_vector_y();
      }
      else {
         // find x-axis which lies on plane with normal vector_z
         // (there must exist a vector with coor {1,0,x} normal to z_vector )
         vdc::vector_x = {1.0, 0.0, -1.0*vdc::vector_z[0]/vdc::vector_z[2]};

         // find vector_y
         init_vector_y();
      }
   }
   else if ( !z_vector && x_vector ){

      // x-axis cannot be initialised alone
      ////terminaltextcolor(RED);
      std::cerr << "Error - x-axis cannot be initialised alone."
                << "\n" << "To use 1D colour scheme, initialise z-axis instead"
                << std::endl;
      ////terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }
   else if ( z_vector && x_vector ){

      // check if input axes are orthogonal
      double zdotx;
      zdotx = vdc::vector_z[0]*vdc::vector_x[0] + vdc::vector_z[1]*vdc::vector_x[1] + vdc::vector_z[2]*vdc::vector_x[2];

      if ( (zdotx > 0.000001) || (zdotx < -0.000001) ){
         ////terminaltextcolor(RED);
         std::cerr << "Error - input axes are not orthogonal." << std::endl;
         ////terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

   }
   return EXIT_SUCCESS;

}

//------------------------------------------------------------------------------
// Extracts 3D vector coordinates from string: x,y,z
// where x,y and z are type double
//------------------------------------------------------------------------------
int extract_vector( std::string arg_string, std::vector<double>& arg_vector ){
   std::string tmp_string;
   int vector_index = 0;

   // read coordinates
   for (unsigned int i = 0; i < arg_string.size(); i++){
      if ( arg_string[i] == ','){
         // check if a number has bean read (no leading comma)
         if ( tmp_string.size() == 0 ){
            std::cerr << "Error - vector should be in the format x,y,z with no spaces." << std::endl;
            return EXIT_FAILURE;
         }

         // save coordinate and move onto next one
         arg_vector[vector_index] = std::stod(tmp_string);
         tmp_string = "";
         vector_index++;
      }
      else if ( i == (arg_string.size()-1) ){
         // reached end of char, read final coordinate
         tmp_string.push_back(arg_string[i]);
         arg_vector[vector_index] = std::stod(tmp_string);
      }
      else{
         tmp_string.push_back(arg_string[i]);
      }
   }

   // check vector has been input correctly (no missing coordinates)
   if ( vector_index != 2){
      std::cerr << "Error - vector should be in the format x,y,z with no spaces." << std::endl;
   }

   // normalise arg_vector
   double length;
   length = std::sqrt( arg_vector[0]*arg_vector[0] + arg_vector[1]*arg_vector[1] + arg_vector[2]*arg_vector[2] );
   arg_vector[0] = arg_vector[0]/length;
   arg_vector[1] = arg_vector[1]/length;
   arg_vector[2] = arg_vector[2]/length;

   return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
// Extracts fractional min max for slice of system to be shown
// xmin,xmax,ymin,ymax,zmin,zmax all type double
//------------------------------------------------------------------------------
int extract_slice_param( std::string arg_string, std::vector<double>& arg_vector, int number_of_param){
   std::string tmp_string = "";
   int vector_index = 0;

   // read fractional coordinates
   for (unsigned int i = 0; i < arg_string.size(); i++){
      if ( arg_string[i] == ',' ){
         // check if a number has bean read (no leading comma)
         if ( tmp_string.size() == 0 ){
            std::cerr << "Error - leading comma." << std::endl;
            return EXIT_FAILURE;
         }

         // save coordinate, verify and move onto next one
         arg_vector[vector_index] = std::stod(tmp_string);

         if ( (arg_vector[vector_index] <= -0.000001) || (arg_vector[vector_index] >= 1.000001 ) ){
            std::cerr << "Error - fractional coordinates should be between 0-1." << std::endl;
            return EXIT_FAILURE;
         }
         tmp_string = "";
         vector_index++;
      }
      else if ( i == (arg_string.size()-1) ){
         // reached end of character, read final number
         tmp_string.push_back(arg_string[i]);
         arg_vector[vector_index] = std::stod(tmp_string);

         if ( (arg_vector[vector_index] <= -0.000001) || (arg_vector[vector_index] >= 1.000001 ) ){
            std::cerr << "Error - fractional coordinates should be between 0-1." << std::endl;
            return EXIT_FAILURE;
         }
      }
      else {
         // push back the digit
         tmp_string.push_back(arg_string[i]);
      }
   }

   // check vector has been input correctly (no missing coordinates)
   if ( vector_index != (number_of_param-1)){
      std::cerr << "Error - 6 comma separated values require for slice or slice-void, 3 for slice-sphere." << std::endl;
   }

   return EXIT_SUCCESS;
}
//------------------------------------------------------------------------------
// Extracts materials to be remove from the visualised system
//------------------------------------------------------------------------------
int extract_materials( std::string arg_string, std::vector<int>& arg_vector){
   std::string tmp_string = "";

   // read materials
   for (unsigned int i = 0; i < arg_string.size(); i++){
      if ( arg_string[i] == ',' ){
         // check if a number has bean read (no leading comma)
         if ( tmp_string.size() == 0 ){
            std::cerr << "Error - leading comma" << std::endl;
            return EXIT_FAILURE;
         }

         // save material and move onto next one
         arg_vector.push_back(std::stod(tmp_string));
         tmp_string = "";
      }
      else if ( i == (arg_string.size()-1) ){
         // reached end of character, read final material
         tmp_string.push_back(arg_string[i]);
         arg_vector.push_back(std::stod(tmp_string));

      }
      else {
         // push back the digit
         tmp_string.push_back(arg_string[i]);
      }
   }

   return EXIT_SUCCESS;
}
//------------------------------------------------------------------------------
// Perform cross product of input vectors vector_x and vector_z to get vector_y
//------------------------------------------------------------------------------
void init_vector_y(){

   vdc::vector_y[0] = vdc::vector_z[1]*vdc::vector_x[2] - vdc::vector_x[1]*vdc::vector_z[2];
   vdc::vector_y[1] = vdc::vector_x[0]*vdc::vector_z[2] - vdc::vector_z[0]*vdc::vector_x[2];
   vdc::vector_y[2] = vdc::vector_z[0]*vdc::vector_x[1] - vdc::vector_x[0]*vdc::vector_z[1];

   return;
}


//------------------------------------------------------------------------------
// Check number of command line args not exceeded
//------------------------------------------------------------------------------
int check_arg( int& arg, int argc, char* argv[], std::string& temp_str, std::string error_output ){

   if (arg+1 < argc){
      arg++;
      temp_str = argv[arg];
   }
   else {
      std::cerr << error_output << std::endl;

      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;

}

} // end of namespace vdc
