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
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations

//------------------------------------------------------------------------------
// Function to output crystal.xyz file compatible with rasmol
//------------------------------------------------------------------------------
void output_xyz_file(){

   // output informative message to user
   if(vdc::verbose) std::cout << "Writing xyz file... " << std::flush;

   // output xyz file
   std::ofstream ofile;
   ofile.open("crystal.xyz");

   // output number of atoms
   ofile << vdc::num_atoms << "\n\n";

   for(uint64_t atom = 0; atom < vdc::num_atoms; atom++){

      // get atom type
      int type_id = vdc::type[atom];

      ofile << materials[type_id].element << "\t" <<
               vdc::coordinates[3*atom + 0] << "\t" <<
               vdc::coordinates[3*atom + 1] << "\t" <<
               vdc::coordinates[3*atom + 2] << "\n";

   }

   ofile << std::flush;
   ofile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
