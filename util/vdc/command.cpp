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
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
int extract( std::string arg_string, std::vector<double>& arg_vector );
int check_arg( int& arg, int argc, char* argv[], std::string& temp_str, std::string error_output );
void init_vector_y(std::vector<double> vector_z, std::vector<double> vector_x );

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
      // Check for colour mapping parameters
      //------------------------------------------------------------------------
      else if (sw == "--vector-z"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables in brackets." );

         // work through vector and extract values
         extract(temp_str, vdc::vector_z );

         // confirm initialisation of z-axis
         z_vector = true;
      }
      else if (sw == "--vector-x"){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables in brackets." );

         // work through vector and extract values
         extract(temp_str, vdc::vector_x );

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
      else if ( sw == "--custom_colourmap" ){

         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - expected custom colourmap name.");

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
         init_vector_y( vdc::vector_z, vdc::vector_x );
      }
      else {
         // find x-axis which lies on plane with normal vector_z
         // (there must exist a vector with coor {1,0,x} normal to z_vector )
         vdc::vector_x = {1.0, 0.0, -1.0*vdc::vector_z[0]/vdc::vector_z[2]};

         // find vector_y
         init_vector_y( vdc::vector_z, vdc::vector_x );
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
// Extracts 3D vector coordinates from string: {x,y,z} or (x,y,z)
// where x,y and z are type double
//------------------------------------------------------------------------------
int extract( std::string arg_string, std::vector<double>& arg_vector ){
   int marker = 0; // position in the vector string

   // check for opening brackets
   if ( (arg_string[marker] != '(') && (arg_string[marker] != '{') ){
      ////terminaltextcolor(RED);
      std::cerr << "Error - brackets required around 3 comma separated values"
                << std::endl;
      ////terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }

   // move to next character
   marker++;

   // read coordinates
   for ( int i = 0; i < 3; i++){
      std::string tmp_string;

      // read coordinate-value
      int j = 0;

      while ( (arg_string[marker] != ',') && (arg_string[marker] != '}') && (arg_string[marker] != ')') ){
         tmp_string.push_back(arg_string[marker]);

         // move through number
         marker++;
         j++;
      }

      arg_vector[i] = std::stod(tmp_string);

      // skip comma, check for closing brackets
      if ( arg_string[marker] == ',' ){
         marker++;
      }
      else if ( ((arg_string[marker] != ')') && (arg_string[marker] != '}' )) && ( i == 2 ) ){
         ////terminaltextcolor(RED);
         std::cerr << "Error - brackets required around 3 comma separated values"
                   << std::endl;
         ////terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
      else if ( ((arg_string[marker] == ')') || (arg_string[marker] == '}')) && ( i != 2 ) ){
         ////terminaltextcolor(RED)
         std::cerr << "Error - three coordinates required"
                   << std::endl;
         ////terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

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
// Perform cross product of input vectors vector_x and vector_z to get vector_y
//------------------------------------------------------------------------------
void init_vector_y(std::vector<double> vector_z, std::vector<double> vector_x ){

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
