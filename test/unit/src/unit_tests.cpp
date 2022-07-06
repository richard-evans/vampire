//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// include header for test functions
#include "unit_tests.hpp"

int main(int argc, char* argv[]){

   // set up module struct with list of tests to be run
   ut::module_t module;

   // specify verbosity of output
   bool verbose = true;

   // process command line options
   //ut::cmd(argc, argv, module, verbose);

   // specify all tests to be run unless arguments are provided
   const bool all = argc > 1 ? false : true;

   // calculate total number of errors
   int error_count = 0;

   //---------------------------------------------------------------------------
   // Run tests
   //---------------------------------------------------------------------------

   std::cout << "--------------------------------------------------" << std::endl;
   std::cout << "    Running unit test suite for vampire code" << std::endl;
   std::cout << "--------------------------------------------------" << std::endl;

   if( module.utility || all ) error_count += ut::utility_tests(verbose);


   // Summary
   std::cout << "--------------------------------------------------" << std::endl;
   if(error_count >0){
      std::cout << "Failed " << error_count << " tests : OVERALL FAIL" << std::endl;
   }
   else{
      std::cout << "Failed " << error_count << " tests : OVERALL PASS" <<	std::endl;
   }
   std::cout << "--------------------------------------------------" << std::endl;

   return 0;

}
