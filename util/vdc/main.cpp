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

// program header
#include "vdc.hpp"

int main(int argc, char* argv[]){

   // process command line arguments
   //command(argc, argv);

   if(vdc::verbose){
      std::cout << "|------------------------------------------------------------|" << std::endl;
      std::cout << "|              Vampire Data Converter for v5+                |" << std::endl;
      std::cout << "|------------------------------------------------------------|" << std::endl;
   }

   // process coordinates
   vdc::process_coordinates();

   vdc::initialise_cells();

   // process spin files
   vdc::process_spins();

   return 0;

}
