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

      }

   }

   void extract( std::string string_z, std::vector<double> vector_z ){
      int marker = 0; //position in the vector string

      // check for opening brackets
      if ( string_z[marker] == ( "(" || "{" ) ) {
         //-----------------------------------------------
         // move to next character
         marker++;

         // read coordinates
         read_vector( vect, x, marker );
         read_vector( vect, y, marker );
         read_vector( vect, z, marker );

         // normalise (x,y,z)
         double length;
         length = std:sqrt( x*x + y*y + z*z );
         x = x/length;
         y = y/length;
         z = z/length;


      }
      else {
         terminaltextcolor(RED);
         std::cerr << "Error - brackets required around 3 comma separated values"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
   }

   void read_vector( std::string vect, double& coordinate, int& marker ){
      std::string tmp_string;
      int i = 0;

      // read coordinate-value
      while ( vect[marker] != "," ) {
         tmp_string.resize( tmp_string.size() +1 );
         tmp_string[i] = vect[marker];
         marker++;
         i++;
      }

      // move marker off comma
      marker++;

      // convert from string to double
      coordinate = std::stod(tmp_string)
   }




} // end of namespace vdc
