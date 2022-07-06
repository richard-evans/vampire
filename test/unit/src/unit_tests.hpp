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

// include header for test functions
#pragma once

namespace ut{

   // simple struct specifying modules to test
   struct module_t {
      bool utility = false;
   };

   // module level functions
   int utility_tests(const bool verbose);

}
