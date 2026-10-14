//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>

// module headers
#include "internal.hpp"

namespace{

// run vampire with the given input file and return time and spin direction from the last output line
bool run_single_spin(const std::string executable, const std::string input_file, double& t, double s[3]){
   if( vt::system(executable + " --input-file " + input_file) != 0 ) return false;
   std::ifstream ifile("output");
   std::string line;
   std::string last_line;
   while(getline(ifile, line)){
      if(!line.empty() && line[0] != '#') last_line = line;
   }
   ifile.close();
   vt::remove({"output", "log"});
   std::stringstream liness(last_line);
   liness >> t >> s[0] >> s[1] >> s[2];
   return !last_line.empty();
}

}

//------------------------------------------------------------------------------
// Test of the order of convergence of an integrator for a single spin with
// damping alpha precessing in a field B along z, starting along x. The exact
// solution at time t is
//
//    S = ( cos(wt)/cosh(alpha wt), sin(wt)/cosh(alpha wt), tanh(alpha wt) ),
//    w = gamma B / (1 + alpha^2)
//
// The data directory holds two inputs to the same final time with time steps
// dt (input-coarse) and dt/2 (input-fine). The measured order
// p = log2(error_coarse/error_fine) must agree with the expected order within tol.
//------------------------------------------------------------------------------
bool convergence_test(const std::string dir, double alpha, double field, double expected_order, double tol, const std::string executable){

   const double gamma = 1.76085963023e11; // gyromagnetic ratio used by vampire (1/(Ts))

   // get root directory
   std::string path = std::filesystem::current_path();

   // fixed-width output for prettiness
   std::stringstream test_name;
   test_name << "Testing convergence order for " << dir;
   std::cout << std::setw(60) << std::left << test_name.str() << " : " << std::flush;

   // change directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;

   double t_coarse = 0.0;
   double t_fine = 0.0;
   double s_coarse[3] = {0.0, 0.0, 0.0};
   double s_fine[3] = {0.0, 0.0, 0.0};
   const bool ran = run_single_spin(executable, "input-coarse", t_coarse, s_coarse) &&
                    run_single_spin(executable, "input-fine",   t_fine,   s_fine);

   // return to parent directory
   if( !vt::chdir(path) ) return false;

   if( !ran || std::abs(t_coarse - t_fine) > 1.0e-6 * t_fine ){
      std::cout << "FAIL | error running vampire or different final times " << t_coarse << " " << t_fine << std::endl;
      return false;
   }

   // exact solution at the final time
   const double wt = gamma * field * t_fine / (1.0 + alpha * alpha);
   const double exact[3] = { std::cos(wt) / std::cosh(alpha * wt), std::sin(wt) / std::cosh(alpha * wt), std::tanh(alpha * wt) };

   double error_coarse = 0.0;
   double error_fine = 0.0;
   for(int i = 0; i < 3; ++i){
      error_coarse += (s_coarse[i] - exact[i]) * (s_coarse[i] - exact[i]);
      error_fine   += (s_fine[i]   - exact[i]) * (s_fine[i]   - exact[i]);
   }
   error_coarse = std::sqrt(error_coarse);
   error_fine   = std::sqrt(error_fine);

   const double order = std::log2(error_coarse / error_fine);

   if( std::abs(order - expected_order) <= tol ){
      std::cout << "OK" << std::endl;
      return true;
   }
   else{
      std::cout << "FAIL | expected order: " << expected_order << "\tobtained: " << order
                << " (errors " << error_coarse << ", " << error_fine << ")" << std::endl;
      return false;
   }

}
