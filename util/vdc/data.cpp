//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// program header
#include "vdc.hpp"

namespace vdc{

   // program option flags
   bool verbose = true; // flag to specify verbosity of output to user

   format_t format;

   uint64_t num_atoms = 0;

   std::vector<material_t> materials(0);

   std::vector<double> coordinates(0);
   std::vector<int> category(0);
   std::vector<int> type(0);
   
   std::vector<double> spins(0);

} // end of namespace vdc
