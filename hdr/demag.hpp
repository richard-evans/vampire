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

#ifndef DEMAG_H_
#define DEMAG_H_

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

namespace demag{

   extern bool fast;
   extern std::uint64_t update_rate;
   extern std::uint64_t update_time;

   extern void init();
   extern void update();

}

#endif //DEMAG_H_
