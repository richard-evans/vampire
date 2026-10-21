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
#include "mpi_tests.hpp"
#include "shmem_test.hpp"

namespace ut{
int mpi_tests(const bool verbose){

   if(verbose) std::cout << "Testing mpi module" << std::endl;

   int error_count = 0;

   error_count += ut::mpi::test_shared_array(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " mpi                  : PASS " << std::endl;
   else std::cout <<                 " mpi                  : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
