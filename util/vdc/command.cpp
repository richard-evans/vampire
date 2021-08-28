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

         // change temp_str to lowercase for comparison
         for (char &c : temp_str){
            c = std::tolower(c);
         }

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
