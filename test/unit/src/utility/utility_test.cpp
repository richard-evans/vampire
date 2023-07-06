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
#include "utility_test.hpp"

namespace ut{
//------------------------------------------------------------------------------
// Function to test utility module functions
//------------------------------------------------------------------------------
int utility_tests(const bool verbose){

   if(verbose) std::cout << "Testing utility module" << std::endl;

   int error_count = 0;

   error_count += ut::utility::test_units(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " utility             : PASS " << std::endl;
   else std::cout <<                 " utility             : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
