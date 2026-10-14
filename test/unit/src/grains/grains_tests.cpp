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
#include "grains_tests.hpp"
#include "grain_geometry_test.hpp"
#include "grain_packing_test.hpp"
#include "grain_seeds_test.hpp"
#include "grain_shape_test.hpp"
#include "grain_tessellation_test.hpp"

namespace ut{
//------------------------------------------------------------------------------
// Function to test grains module functions
//------------------------------------------------------------------------------
int grains_tests(const bool verbose){

   if(verbose) std::cout << "Testing grains module" << std::endl;

   int error_count = 0;

   error_count += ut::grains::test_grain_geometry(verbose);
   error_count += ut::grains::test_grain_tessellation(verbose);
   error_count += ut::grains::test_grain_packing(verbose);
   error_count += ut::grains::test_grain_seeds(verbose);
   error_count += ut::grains::test_grain_shape(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " grains               : PASS " << std::endl;
   else std::cout <<                 " grains               : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
