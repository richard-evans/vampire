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
void extract( std::string string_z, std::vector<double> vector_z );

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
   std::string string_z = "(0.0,0.0,1.0)"
   std::string string_x = "(1.0,0.0,0.0)"

   std::vector<double> vector_z = {0.0,0.0,1.0};
   std::vector<double> vector_x = {1.0,0.0,0.0};

   for (int arg = 1; arg < argc; arg++){

      // read prefix
      std::string sw=argv[arg]
      if (sw == "--vector-z"){

         // check number of args not exceeded
         if (arg+1 < argc){
            arg++;
            vector_z = string(argv[arg]);
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - expected 3 comma separated variables in brackets."
                      << "\n" << "Check for spaces in command-line arguments"
                      << std:endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;
         }

         // work through vector and extract values
         extract(string_z,vector_z);
      }
      else if (sw == "--vector-x"){

         // check number of args not exceeded
         if (arg+1 < argc){
            arg++;
            vector_z = string(argv[arg]);
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - expected 3 comma separated variables in brackets."
                      << "\n" << "Check for spaces in command-line arguments"
                      << std:endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;
         }

         // work through vector and extract values
         extract(string_z,vector_x);
      }

   }

}

//------------------------------------------------------------------------------
// Extracts 3D vector coordinates from string: {x,y,z} or (x,y,z)
// where x,y and z are doubles
//------------------------------------------------------------------------------
void extract( std::string string_z, std::vector<double> vector_z ){
   int marker = 0; // position in the vector string

   // check for opening brackets
   if ( string_z[marker] != ( "(" || "{" ) ){
      terminaltextcolor(RED);
      std::cerr << "Error - brackets required around 3 comma separated values"
                << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }

   // move to next character
   marker++;

   // read coordinates
   for ( int i = 0; i < 3; i++){
      std::string tmp_string;

      // read coordinate-value
      int j = 0;
      while ( string_z[marker] != "," ){
         tmp_string.resize( tmp_string.size() +1 );
         tmp_string[j] = string_z[marker];

         // move through number
         marker++;
         j++;
      }

      vector_z[i] = std::stod(tmp_string);

      // skip comma, check for closing brackets
      if ( string_z[marker] == "," ){
         marker++;
      }
      else if ( (string_z[marker] != ( ")" || "}" )) && ( i == 2 ){
         terminaltextcolor(RED);
         std::cerr << "Error - brackets required around 3 comma separated values"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
      else if ( (string_z[marker] == ( ")" || "}" )) && ( i != 2 ){
         terminaltextcolor(RED)
         std::cerr << "Error - three coordinates required"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

   }

   // normalise vector_z
   double length;
   length = std:sqrt( vector_z[0]*vector_z[0] + vector_z[1]*vector_z[1] + vector_z[2]*vector_z[2] );
   vector_z[0] = vector_z[0]/length;
   vector_z[1] = vector_z[1]/length;
   vector_z[2] = vector_z[2]/length;

}

} // end of namespace vdc
