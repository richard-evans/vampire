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
#include <cmath>

// module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
// Test to verify consistent time evolution for different integrators
//     Note: this is a regression test rather than absolute as some tests
//           are stochastic in nature
//
//     The spin direction is read from the output line at the given time
//     (first column) and each component must agree with the expected value
//     within |obtained - expected| <= abs_tol + rel_tol*|expected|
//------------------------------------------------------------------------------
bool integrator_test(const std::string dir, double time, double rx, double ry, double rz, const std::string executable){

   // tolerances, set by the 6 significant figures written to the output file
   const double rel_tol = 1.0e-5;
   const double abs_tol = 1.0e-6;

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
      vt::remove({"output", "log"});
      vt::chdir(path);
      return false;
   }

   // open output file
   std::ifstream ifile;
   ifile.open("output");

   // find the data line at the requested time, skipping header lines
   std::string line;
   std::string time_line;
   while(getline(ifile, line)){
      if(line.empty() || line[0] == '#') continue;
      std::stringstream liness(line);
      double t = 0.0;
      liness >> t;
      if( std::abs(t - time) <= 1.0e-6 * std::abs(time) ){
         time_line = line;
         break;
      }
   }
   ifile.close();

   double vx = 0.0;
   double vy = 0.0;
   double vz = 0.0;
   double t = 0.0;
   std::stringstream liness(time_line);
   liness >> t >> vx >> vy >> vz;

   // cleanup
   vt::remove({"output", "log"});

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   if( time_line.empty() ){
      std::cout << "FAIL | no output line at time " << time << std::endl;
      return false;
   }

   // now test value obtained from code
   const bool okx = std::abs(vx - rx) <= abs_tol + rel_tol * std::abs(rx);
   const bool oky = std::abs(vy - ry) <= abs_tol + rel_tol * std::abs(ry);
   const bool okz = std::abs(vz - rz) <= abs_tol + rel_tol * std::abs(rz);

   if( okx && oky && okz ){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected: " << rx << "\t" << ry << "\t" << rz << "\t" << "\tobtained:  " << vx << "\t" << vy << "\t" << vz << "\t" << "\t" << time_line << std::endl;
      return false;
   }

}
