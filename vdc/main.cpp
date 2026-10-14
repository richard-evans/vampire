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
#include <string>

// program header
#include "vdc.hpp"

int main(int argc, char* argv[]){

   std::cout << "|------------------------------------------------------------|" << std::endl;
   std::cout << "|              Vampire Data Converter for v5+                |" << std::endl;
   std::cout << "|------------------------------------------------------------|" << std::endl;

   // process command line arguments
   vdc::command(argc, argv);

   // process input file
   vdc::read_and_set();

   // process coordinates
   vdc::process_coordinates();

   // process spin files
   vdc::process_spins();

   return EXIT_SUCCESS;
}
