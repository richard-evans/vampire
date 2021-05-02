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
#include <algorithm>
#include <fstream>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
// void extract_vector( std::string arg_string, std::vector<double>& arg_vector );
// void extract_materials( std::string arg_string, std::vector<int>& arg_vector);
// void extract_slice_param( std::string arg_string, std::vector<double>& arg_vector, int number_of_param);
void check_arg( int& arg, int argc, char* argv[], std::string& temp_str, std::string error_output );
void init_vector_y();

//------------------------------------------------------------------------------
// Command line parsing function
//------------------------------------------------------------------------------
void command( int argc, char* argv[] ){
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
      if      (sw == "--xyz"   ){ vdc::xyz    = true; } // xyz coordinate file output
      else if (sw == "--povray"){ vdc::povray = true; } // pov coordinate file output
      else if (sw == "--vtk"   ){ vdc::vtk    = true; } // vtk coordinate file output
      else if (sw == "--text"  ){ vdc::txt    = true; } // plain text file output
      else if (sw == "--cells" ){ vdc::cells  = true; }
      else if (sw == "--ssc" || sw == "--spin-spin-correlation"){ vdc::ssc = true; }
      //------------------------------------------------------------------------
      // Check for cell size
      //------------------------------------------------------------------------
      else if (sw == "--cell-size"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - size of cells in Angstroms." );

         // set cell size
         if ( stof(temp_str) >= 0.5 ){ vdc::cell_size = stof(temp_str); }
         else{
            std::cerr << "Error - cell size must be greater than 0.5 Angstroms." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }
      //------------------------------------------------------------------------
      // Check for verbose output
      //------------------------------------------------------------------------
      else if (sw == "--verbose" || sw == "-v"){ vdc::verbose = true; }
      //------------------------------------------------------------------------
      // Input file name
      //------------------------------------------------------------------------
      else if (sw == "--input-file"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - no file specified for \'--input-file\' command line option." );

         // set input_file name
         vdc::input_file = temp_str;
      }
      //------------------------------------------------------------------------
      // Help information
      //------------------------------------------------------------------------
      else if (sw == "--help" || sw == "-h"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - no parameter specified for "+sw+" command line option." );

         // check if parameter exists
         if (vdc::key_list.find(temp_str) == vdc::key_list.end()){
            std::cerr << "Error - No matching parameter '" << temp_str << "'.\n";
            std::exit(EXIT_FAILURE);
         }

         input_t input;

         // create argument to pass to function wrapper
         input.value = {"-h"};

         // pass to function wrapper
         vdc::key_list.at(temp_str)(input);
      }
      // else if (sw == "--gen-input" || sw == "--generate-input"){

      //    // open input file to write to
      //    std::ofstream input_file(vdc::input_file);

      //    // check if file open was success
      //    if (!input_file.is_open()){
      //       std::cerr << "Error - Unable to open '" << vdc::input_file << "'.\n";
      //       std::exit(EXIT_FAILURE);
      //    }

      //    input_file << "#=================================\n"
      //               << "# Povray Parameters\n"
      //               << "#=================================\n\n"
      //               << "colourmap = CBWR\n"
      //               << ""
                    



      //    std::exit(EXIT_SUCCESS);
      // }
      //------------------------------------------------------------------------
      // Check for user specified initial frame to render
      //------------------------------------------------------------------------
      // else if (sw == "--frame-start"){
      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected index of intial frame to render." );

      //    if ( stoi(temp_str) >= 0 ){ vdc::vdc_start_file_id = stoi(temp_str); }
      //    else {
      //       std::cerr << "Error - frame index cannot be negative." << std::endl;
      //       std::exit(EXIT_FAILURE);
      //    }
      // }
      //------------------------------------------------------------------------
      // Check for user specified final frame to render
      //------------------------------------------------------------------------
      // else if (sw == "--frame-final"){
      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected index of final frame to render." );

      //    // +1 since we want to render the last frame as well
      //    if ( stoi(temp_str) >= 0 ){ vdc::vdc_final_file_id = stoi(temp_str)+1; }
      //    else {
      //       std::cerr << "Error - frame index cannot be negative." << std::endl;
      //       std::exit(EXIT_FAILURE);
      //    }
      // }
      //------------------------------------------------------------------------
      // Check for user specified materials to remove
      //------------------------------------------------------------------------
      // else if (sw == "--remove-material" || sw == "--remove-materials"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected at least one variable." );

      //    // work through vector and extract values
      //    extract_materials(temp_str, vdc::remove_materials);
      // }
      //------------------------------------------------------------------------
      // Check for user specified antiferromagnetic material spins to flip in
      // povray visualisation
      //------------------------------------------------------------------------
      // else if (sw == "--afm"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected at least one variable." );

      //    // work through vector and extract values
      //    extract_materials(temp_str, vdc::afm_materials);
      // }
      //------------------------------------------------------------------------
      // Check for slice parameters
      //------------------------------------------------------------------------
      // else if (sw == "--slice"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 6 comma separated variables." );

      //    // work through vector and extract values
      //    extract_slice_param(temp_str, vdc::slice_parameters, 6);

      //    // set slice keyword
      //    vdc::slice_type = vdc::slice;
      // }
      // else if (sw == "--slice-void"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 6 comma separated variables." );

      //    // work through vector and extract values
      //    extract_slice_param(temp_str, vdc::slice_parameters, 6);

      //    // set slice keyword
      //    vdc::slice_type = vdc::slice_void;
      // }
      // else if (sw == "--slice-sphere"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

      //    // work through vector and extract values
      //    extract_slice_param(temp_str, vdc::slice_parameters, 3);

      //    // set slice keyword
      //    vdc::slice_type = vdc::slice_sphere;
      // }
      // else if (sw == "--slice-cylinder"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 4 comma separated variables." );

      //    // work through vector and extract values
      //    extract_slice_param(temp_str, vdc::slice_parameters, 4);

      //    // set slice keyword
      //    vdc::slice_type = vdc::slice_cylinder;
      // }
      //------------------------------------------------------------------------
      // Check for colour mapping parameters
      //------------------------------------------------------------------------
      // else if (sw == "--vector-z"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

      //    // work through vector and extract values
      //    extract_vector(temp_str, vdc::vector_z );

      //    // confirm initialisation of z-axis
      //    z_vector = true;
      // }
      // else if (sw == "--vector-x"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected 3 comma separated variables." );

      //    // work through vector and extract values
      //    extract_vector(temp_str, vdc::vector_x );

      //    // confirm initialisation of x-axis
      //    x_vector = true;
      // }
      // else if (sw == "--colourmap"){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected colourmap keyword." );

      //    if (std::find(vdc::colourmaps.begin(),vdc::colourmaps.end(),temp_str) != vdc::colourmaps.end() ){
      //       vdc::colour_keyword = temp_str;
      //    }
      //    else {
      //       std::cerr << "Error - Colourmap keyword does not match." << std::endl;
      //       std::exit(EXIT_FAILURE);
      //    }
      // }
      // else if ( sw == "--3D" ){ vdc::x_axis_colour = true; }
      // else if ( sw == "--custom-colourmap" ){

      //    // check number of args not exceeded
      //    check_arg(arg, argc, argv, temp_str, "Error - expected custom colourmap name.");

      //    // set colour_keyword to "custom"
      //    vdc::colour_keyword = "custom";

      //    // set custom map file name
      //    vdc::custom_colourmap_file = temp_str;

      // }
      else {
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
         std::exit(EXIT_FAILURE);
      }
   }

   //---------------------------------------------------------------------------
   // Additional checks on command line parameters
   //---------------------------------------------------------------------------

   // check that some kind of data output is requested
   if( !vdc::xyz && !vdc::povray && !vdc::vtk && !vdc::txt && !vdc::ssc && !vdc::cells){
      std::cerr << "Error! No output data formats requested. Available options are: " << std::endl;
      std::cerr << "\t\t --xyz    Data output in .xyz format for viewing in rasmol/jmol" << std::endl;
      std::cerr << "\t\t --povray Data output in PoVRAY format for rendering" << std::endl;
      std::cerr << "\t\t --vtk    Data output in VTK format for viewing in Paraview" << std::endl;
      std::cerr << "\t\t --text   Data output in plain text format for plotting in gnuplot/excel etc" << std::endl;
      std::cerr << "\t\t --cells  Data output in plain text format in cells" << std::endl;
      std::cerr << "\t\t --ssc    Spin-spin correlation data in text format" << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

//------------------------------------------------------------------------------
// Extracts 3D vector coordinates from string: x,y,z
// where x,y and z are type double
//------------------------------------------------------------------------------
// void extract_vector( std::string arg_string, std::vector<double>& arg_vector ){
//    std::string tmp_string;
//    int vector_index = 0;

//    // read coordinates
//    for (unsigned int i = 0; i < arg_string.size(); i++){
//       if ( arg_string[i] == ','){
//          // check if a number has bean read (no leading comma)
//          if ( tmp_string.size() == 0 ){
//             std::cerr << "Error - vector should be in the format x,y,z with no spaces." << std::endl;
//             std::exit(EXIT_FAILURE);
//          }

//          // save coordinate and move onto next one
//          arg_vector[vector_index] = std::stod(tmp_string);
//          tmp_string = "";
//          vector_index++;
//       }
//       else if ( i == (arg_string.size()-1) ){
//          // reached end of char, read final coordinate
//          tmp_string.push_back(arg_string[i]);
//          arg_vector[vector_index] = std::stod(tmp_string);
//       }
//       else{
//          tmp_string.push_back(arg_string[i]);
//       }
//    }

//    // check vector has been input correctly (no missing coordinates)
//    if ( vector_index != 2){
//       std::cerr << "Error - vector should be in the format x,y,z with no spaces." << std::endl;
//    }

//    // normalise arg_vector
//    double length;
//    length = std::sqrt( arg_vector[0]*arg_vector[0] + arg_vector[1]*arg_vector[1] + arg_vector[2]*arg_vector[2] );
//    arg_vector[0] = arg_vector[0]/length;
//    arg_vector[1] = arg_vector[1]/length;
//    arg_vector[2] = arg_vector[2]/length;
// }

// //------------------------------------------------------------------------------
// // Extracts fractional min max for slice of system to be shown
// // xmin,xmax,ymin,ymax,zmin,zmax all type double
// //------------------------------------------------------------------------------
// void extract_slice_param( std::string arg_string, std::vector<double>& arg_vector, int number_of_param){
//    std::string tmp_string = "";
//    int vector_index = 0;

//    // read fractional coordinates
//    for (unsigned int i = 0; i < arg_string.size(); i++){
//       if ( arg_string[i] == ',' ){
//          // check if a number has bean read (no leading comma)
//          if ( tmp_string.size() == 0 ){
//             std::cerr << "Error - leading comma." << std::endl;
//             std::exit(EXIT_FAILURE);
//          }

//          // save coordinate, verify and move onto next one
//          arg_vector[vector_index] = std::stod(tmp_string);

//          if ( (arg_vector[vector_index] <= -0.000001) || (arg_vector[vector_index] >= 1.000001 ) ){
//             std::cerr << "Error - fractional coordinates should be between 0-1." << std::endl;
//             std::exit(EXIT_FAILURE);
//          }
//          tmp_string = "";
//          vector_index++;
//       }
//       else if ( i == (arg_string.size()-1) ){
//          // reached end of character, read final number
//          tmp_string.push_back(arg_string[i]);
//          arg_vector[vector_index] = std::stod(tmp_string);

//          if ( (arg_vector[vector_index] <= -0.000001) || (arg_vector[vector_index] >= 1.000001 ) ){
//             std::cerr << "Error - fractional coordinates should be between 0-1." << std::endl;
//             std::exit(EXIT_FAILURE);
//          }
//       }
//       else {
//          // push back the digit
//          tmp_string.push_back(arg_string[i]);
//       }
//    }

//    // check vector has been input correctly (no missing coordinates)
//    if ( vector_index != (number_of_param-1)){
//       std::cerr << "Error - 6 comma separated values require for slice or slice-void, 3 for slice-sphere." << std::endl;
//    }
// }
// //------------------------------------------------------------------------------
// // Extracts materials to be remove from the visualised system
// //------------------------------------------------------------------------------
// void extract_materials( std::string arg_string, std::vector<int>& arg_vector){
//    std::string tmp_string = "";

//    // read materials
//    for (unsigned int i = 0; i < arg_string.size(); i++){
//       if ( arg_string[i] == ',' ){
//          // check if a number has bean read (no leading comma)
//          if ( tmp_string.size() == 0 ){
//             std::cerr << "Error - leading comma" << std::endl;
//             std::exit(EXIT_FAILURE);
//          }

//          // save material and move onto next one
//          arg_vector.push_back(std::stod(tmp_string));
//          tmp_string = "";
//       }
//       else if ( i == (arg_string.size()-1) ){
//          // reached end of character, read final material
//          tmp_string.push_back(arg_string[i]);
//          arg_vector.push_back(std::stod(tmp_string));

//       }
//       else {
//          // push back the digit
//          tmp_string.push_back(arg_string[i]);
//       }
//    }
// }

//------------------------------------------------------------------------------
// Check number of command line args not exceeded
//------------------------------------------------------------------------------
void check_arg( int& arg, int argc, char* argv[], std::string& temp_str, std::string error_output ){

   if (arg+1 < argc){
      arg++;
      temp_str = argv[arg];
   }
   else {
      std::cerr << error_output << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

} // end of namespace vdc
