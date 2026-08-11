//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// include header for test functions
#include "create_tests.hpp"
#include "voronoi_bimodal_seeds_test.hpp"

namespace ut{
//------------------------------------------------------------------------------
// Function to test create module functions
//------------------------------------------------------------------------------
int create_tests(const bool verbose){

   if(verbose) std::cout << "Testing create module" << std::endl;

   int error_count = 0;

   error_count += ut::create::test_voronoi_bimodal_seeds(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " create               : PASS " << std::endl;
   else std::cout <<                 " create               : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
