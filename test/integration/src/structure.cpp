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

// module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
// Test to verify consistent generation of different systems like particles,
// core-shell systems, multilayers etc, tested against number of atoms for each
// material
//------------------------------------------------------------------------------
bool material_atoms_test(const std::string dir, int n1, int n2, int n3, int n4, const std::string executable){

   // get root directory
   std::string path = std::filesystem::current_path();

   // fixed-width output for prettiness
   std::stringstream test_name;
   test_name << "Testing structure for " << dir;
   std::cout << std::setw(60) << std::left << test_name.str() << " : " << std::flush;

   // change directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;

   // run vampire
   int vmp = vt::system(executable);
   if( vmp != 0){
      std::cerr << "Error running vampire. Returning as failed test." << std::endl;
      return false;
   }

   std::string line;

   // create temporary file for material numbers
   std::string cmd = "grep Material log | awk '{print $14}' > nums.txt";
   system(cmd.c_str());

   // open output file
   std::ifstream ifile;
   ifile.open("nums.txt");

   int v1 = 0;
   int v2 = 0;
   int v3 = 0;
   int v4 = 0;

   getline(ifile, line);
   std::stringstream liness1(line);
   liness1 >> v1;

   getline(ifile, line);
   std::stringstream liness2(line);
   liness2 >> v2;

   getline(ifile, line);
   std::stringstream liness3(line);
   liness3 >> v3;

   getline(ifile, line);
   std::stringstream liness4(line);
   liness4 >> v4;

   // cleanup
   vt::system("rm output log nums.txt");

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   // now test value obtained from code

   if( n1 == v1 && n2 == v2 && n3 == v3 && n4 == v4){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << n1 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << "\tobtained:  " << v1 << "\t" << v2 << "\t" << v3 << "\t" << v4 << "\t" << line << std::endl;
      return false;
   }

}
