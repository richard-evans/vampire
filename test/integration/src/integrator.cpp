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
// Test to verify consistent time evolution for different integrators
//     Note: this is a regression test rather than absolute as some tests
//           are stochastic in nature
//------------------------------------------------------------------------------
bool integrator_test(const std::string dir, double rx, double ry, double rz, const std::string executable){

   // get root directory
   std::string path = std::filesystem::current_path();

   // fixed-width output for prettiness
   std::stringstream test_name;
   test_name << "Testing integrator for " << dir;
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

   // open output file
   std::ifstream ifile;
   ifile.open("output");

   // read value after header
   for(int i=0; i<982; i++) getline(ifile, line);
   std::stringstream liness(line);
   double v1 = 0.0;
   double vx = 0.0;
   double vy = 0.0;
   double vz = 0.0;
   double vm = 0.0;

   liness >> v1 >> vx >> vy >> vz >> vm;

   // cleanup
   vt::system("rm output log");

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   // now test value obtained from code
   const double ratiox = vx/rx;
   const double ratioy = vy/ry;
   const double ratioz = vz/rz;

   if(ratiox >0.99999 && ratiox < 1.00001 && ratioy >0.99999 && ratioy < 1.00001 && ratioz >0.99999 && ratioz < 1.00001){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << rx << "\t" << ry << "\t" << rz << "\t" << "\tobtained:  " << vx << "\t" << vy << "\t" << vz << "\t" << "\t" << line << std::endl;
      return false;
   }

}
