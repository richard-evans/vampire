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

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

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
   ofile << vdc::sliced_atoms_list.size() + vdc::sliced_nm_atoms_list.size() << "\n\n";

   #pragma omp parallel
   {

      std::stringstream otext;

      // write magnetic atoms to output text stream in parallel
      #pragma omp for
      for(int i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

         // get atom type
         int type_id = vdc::type[atom];

         otext << materials[type_id].element << "\t" <<
                  vdc::coordinates[3*atom + 0] << "\t" <<
                  vdc::coordinates[3*atom + 1] << "\t" <<
                  vdc::coordinates[3*atom + 2] << "\n";

      } // end of parallel for

      // write non-magnetic atoms
      #pragma omp for
      for(int i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_nm_atoms_list[i];

         // get atom type
         int type_id = vdc::nm_type[atom];

         otext << materials[type_id].element << "\t" <<
                  vdc::nm_coordinates[3*atom + 0] << "\t" <<
                  vdc::nm_coordinates[3*atom + 1] << "\t" <<
                  vdc::nm_coordinates[3*atom + 2] << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      ofile << otext.str();

   } // end of parallel region

   ofile << std::flush;
   ofile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
