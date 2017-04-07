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

#ifndef VDC_H_
#define VDC_H_

// C++ standard library headers
#include <iostream>
#include <vector>

namespace vdc{

   // program option flags
   extern bool verbose;
   extern bool xyz;

   // enumerated integers for option selection
   enum format_t{ binary = 0, text = 1};
   extern format_t format;

   // simple struct to store material parameters
   struct material_t{
      int id;
      double moment;
      std::string name;
      std::string element;
   };

   extern uint64_t num_atoms;

   extern std::vector<material_t> materials;

   extern std::vector<double> coordinates;
   extern std::vector<int> category;
   extern std::vector<int> type;

   extern std::vector<double> spins;

   // array to store subsidiary data file names
   extern std::vector <std::string> coord_filenames;

   // Functions
   void process_coordinates();
   void output_xyz_file();

}

#endif //VDC_H_
