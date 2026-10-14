//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2009-2012. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
//====================================================================================================
//
//                                           Demag.h
//
//                         Namespace header file for demag routines
//
//                               Version 1.0 R Evans 18/09/2009
//
//====================================================================================================
#include <vector>
#include <cstdint>

// Load standard forms of fixed-width types (needed for some compilers)
using std::uint64_t;
using std::int64_t;

namespace demag{

   extern bool fast;
   extern uint64_t update_rate;
   extern uint64_t update_time;

   extern void init();
   extern void update();

}
