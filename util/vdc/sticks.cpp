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
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
bool read_spin_metadata(unsigned int file_id);
void read_spin_data();

//------------------------------------------------------------------------------
// Function to output sticks.inc file compatible with povray
//------------------------------------------------------------------------------
void output_sticks_file(){

   // output informative message to user
   if(vdc::verbose) std::cout << "   Generating stick data and writing sticks file sticks.inc..." << std::flush;

   // open incfile
   std::ofstream sfile;
   sfile.open("sticks.inc");

   const double scr2 = vdc::sticks_cutoff*vdc::sticks_cutoff;

   std::stringstream otext; // output file text asembled in memory for speed

   // loop over all atoms to find nearby atoms within cutoff
   for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

      // get atom ID
      unsigned int atomi = vdc::sliced_atoms_list[i];

      // only for i > j
      for(size_t j=i+1; j < vdc::sliced_atoms_list.size(); j++){

         // get atom ID
         unsigned int atomj = vdc::sliced_atoms_list[j];

         // calculate range between atoms j-i
         const double dx = (coordinates[3*atomj+0]-coordinates[3*atomi+0]);
         const double dy = (coordinates[3*atomj+1]-coordinates[3*atomi+1]);
         const double dz = (coordinates[3*atomj+2]-coordinates[3*atomi+2]);
         const double dr2 = dx*dx + dy*dy + dz*dz;

         // check atoms are within cutoff
         if(dr2 < scr2){

            const int mati = vdc::type[atomi]+1;
            const int matj = vdc::type[atomj]+1;

            // stick(sx,sy,sz,ex,ey,ez,ri,rj)
            // calculate start and end of cylinder based on radius

            const double sx = coordinates[3*atomi+0]-vdc::system_centre[0];
            const double sy = coordinates[3*atomi+1]-vdc::system_centre[1];
            const double sz = coordinates[3*atomi+2]-vdc::system_centre[2];

            const double ex = coordinates[3*atomj+0]-vdc::system_centre[0];
            const double ey = coordinates[3*atomj+1]-vdc::system_centre[1];
            const double ez = coordinates[3*atomj+2]-vdc::system_centre[2];

            otext << "stick(" << sx << ", " << sy << ", " << sz << ", " <<
                                 ex << ", " << ey << ", " << ez << ", rscale" << mati << ", rscale" << matj << ")\n";
         }

      }

   }

   // output string to file
   sfile << otext.str();

   // flush data to include file and close
   sfile << std::flush;
   sfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
