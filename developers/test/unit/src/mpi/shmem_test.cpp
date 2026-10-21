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
//   Unit tests for the shared memory array (src/mpi/shmem.cpp) in the serial
//   build, where the array is ordinary memory held by the calling process
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// include header for test functions
#include "vmpi.hpp"
#include "shmem_test.hpp"

namespace ut{
   namespace mpi{

int test_shared_array(const bool verbose){

   int error_count = 0;

   // a serial build has one rank in its region
   vmpi::initialise_shared_memory();
   if(vmpi::my_shared_rank != 0 || vmpi::num_shared_ranks != 1){
      std::cout << "FAIL: shared memory region should hold one rank in serial" << std::endl;
      error_count++;
   }

   // allocate, check zero initialisation, write and read back
   const uint64_t n = 1000;
   vmpi::shared_array_t<double> array;
   array.allocate(n);
   if(array.data() == nullptr || array.size() != n || !array.is_writer()){
      std::cout << "FAIL: shared array not allocated with " << n << " elements" << std::endl;
      return error_count + 1;
   }
   int nonzero = 0;
   for(uint64_t i = 0; i < n; i++) if(array[i] != 0.0) nonzero++;
   if(nonzero > 0){
      std::cout << "FAIL: " << nonzero << " shared array elements not set to zero" << std::endl;
      error_count++;
   }
   for(uint64_t i = 0; i < n; i++) array[i] = double(i);
   array.synchronise();
   int wrong = 0;
   for(uint64_t i = 0; i < n; i++) if(array[i] != double(i)) wrong++;
   if(wrong > 0){
      std::cout << "FAIL: " << wrong << " shared array elements changed by synchronisation" << std::endl;
      error_count++;
   }

   // reallocation replaces the contents
   array.allocate(10);
   if(array.size() != 10 || array[9] != 0.0){
      std::cout << "FAIL: shared array not reset when reallocated" << std::endl;
      error_count++;
   }

   // a second array of another type shares the same communicators
   vmpi::shared_array_t<int> indices;
   indices.allocate(5);
   indices[4] = 7;
   indices.synchronise();
   if(indices[4] != 7){
      std::cout << "FAIL: integer shared array not written" << std::endl;
      error_count++;
   }

   // free resets the array, and an empty array can be allocated
   array.free();
   indices.free();
   if(array.data() != nullptr || array.size() != 0){
      std::cout << "FAIL: shared array not reset when freed" << std::endl;
      error_count++;
   }
   array.allocate(0);
   array.free();
   vmpi::finalise_shared_memory();

   if(verbose) std::cout << "   shared array                 : " << (error_count == 0 ? "PASS" : "FAIL") << std::endl;

   return error_count;

}

}
}
