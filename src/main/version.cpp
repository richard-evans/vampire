//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>

// Vampire headers
#include "info.hpp"

namespace vinfo{

   // variable string to store version number of code
   std::string vampire_version = "6.0.0"; // vampire code version

   // wrapper function to return vampire version
   std::string version(){
   	return vampire_version;
   }

} // end of vinfo namespace
