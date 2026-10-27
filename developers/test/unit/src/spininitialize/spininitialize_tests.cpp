//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// unit test headers
#include "spininitialize_tests.hpp"
#include "domain_wall_test.hpp"

namespace ut{
int spininitialize_tests(const bool verbose){

   if(verbose) std::cout << "Testing spininitialize module" << std::endl;

   int error_count = 0;

   error_count += ut::spininitialize::test_domain_wall(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " spininitialize       : PASS " << std::endl;
   else std::cout <<                 " spininitialize       : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
