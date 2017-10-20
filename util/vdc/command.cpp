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
      std::string vector_z = "(0,0,1)"
      std::string vector_x = "(1,0,0)"

      double x1,y1,z1,x2,y2,z2;

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
            extract(vector_z,x1,y1,z1);
         }

      }

   }

   void extract( std::string vect, double& x, double& y, double& z );

      // check for opening brackets
      if ( vect[0] == ( "(" || "{" ) ) {
         int i = 1;

         std::string tmp_string;
         // read x-value
         while ( vect[i] != "," ) {
            tmp_string.resize( tmp_string.size() +1 );
            tmp_string[i-1] = vect[i];
         }

         // convert from string to double
         x = std::stod(tmp_string)
      }
      else {
         terminaltextcolor(RED);
         std::cerr << "Error - brackets required around 3 comma separated values"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

      

} // end of namespace vdc
