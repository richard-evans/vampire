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
#include "integrate_tests.hpp"
#include "rotations_test.hpp"

namespace ut{
int integrate_tests(const bool verbose){

   if(verbose) std::cout << "Testing integrate module" << std::endl;

   int error_count = 0;

   error_count += ut::integrate::test_rotations(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " integrate            : PASS " << std::endl;
   else std::cout <<                 " integrate            : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
