//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Function to output spins-xxxxxxxx.txt file in plaintext format
//------------------------------------------------------------------------------
//
// Writes a single file for each snapshot ID formatted as:
//  cx   cy   cz   sx   sy   sz
//
//------------------------------------------------------------------------------
void output_txt_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream txt_file_sstr;
	txt_file_sstr << "spins-";
	txt_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	txt_file_sstr << ".txt";
	std::string txt_file = txt_file_sstr.str();

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing text file " << txt_file << "..." << std::flush;

   // open incfile
   std::ofstream txtfile;
   txtfile.open(txt_file.c_str());

   // output number of atoms
   txtfile << vdc::sliced_atoms_list.size() + vdc::sliced_nm_atoms_list.size() << "\n";

   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   // step 1: parallel formatted output to stringstream in memory
   // step 2: binary write of formatted text to output file (awesomely fast!)
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

         // format text for plain text file
         otext << coordinates[3*atom+0]-vdc::system_centre[0] << "\t" << coordinates[3*atom+1]-vdc::system_centre[1] << "\t" << coordinates[3*atom+2]-vdc::system_centre[2] << "\t" <<
                  spins[3*atom+0] << "\t" << spins[3*atom+1] << "\t" << spins[3*atom+2] << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      txtfile << otext.str();

   } // end of parallel region

   //---------------------------------------------------------------------------
   // write non-magnetic atoms to txt file
   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for(size_t i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_nm_atoms_list[i];

         // format text for text file
         otext << nm_coordinates[3*atom+0]-vdc::system_centre[0] << "\t" << nm_coordinates[3*atom+1]-vdc::system_centre[1] << "\t" << nm_coordinates[3*atom+2]-vdc::system_centre[2] << "\t" <<
                  0.0 << "\t" << 0.0 << "\t" << 0.0 << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      txtfile << otext.str();

   } // end of parallel region

   // flush data to include file and close
   txtfile << std::flush;
   txtfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
