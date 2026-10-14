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
// Overview
// --------
// This file contains integration tests for the Monte Carlo integrator.
//
// The key difference from the LLG integrator tests (integrator.cpp) is that
// Monte Carlo dynamics are stochastic even at T=0 (each trial move is drawn
// from a random distribution), so the exact spin trajectory is not
// reproducible across platforms or compiler versions. Instead, these tests
// verify physical *outcomes* that must hold for any correct implementation:
// for example, that a ferromagnet reaches full saturation at T=0 regardless
// of the starting configuration.
//
// How the test works
// ------------------
// final_value_test() runs vampire in a subdirectory, reads the LAST data
// line of the output file, extracts a single scalar, and checks it matches
// the expected value to within a relative tolerance.
//
// Reading the last line (rather than a fixed line number) means the test is
// insensitive to changes in the number of output steps, making it robust to
// future changes in the input file.
//
// Adding new Monte Carlo tests
// ----------------------------
// 1. Create a subdirectory under developers/test/integration/data/ with an input file
//    and material file.
// 2. Choose a physical observable whose final value has a known analytic
//    result (e.g. magnetisation length = 1 for a T=0 ferromagnet).
// 3. Call final_value_test() from main.cpp with the directory, expected
//    value, and tolerance. To check several output columns from a single
//    simulation, use final_values_test() with a list of expected values.
//
//------------------------------------------------------------------------------

// module headers
#include "internal.hpp"

//------------------------------------------------------------------------------
// final_values_test
//
// Runs a vampire simulation in the given subdirectory, reads the last data
// line of the output file, and checks one or more columns against expected
// values, each to within its own tolerance
// (i.e. |value - expected| <= abs_tol + tol*|expected|). Each value is
// reported on its own line.
//
// Parameters:
//   dir        - path relative to developers/test/integration/data/
//   values     - list of {label, column, expected, tol[, abs_tol]} to check,
//                where column is the 0-based index of the column in the output
//                and abs_tol (default 0) is needed when expected is zero
//   executable - full shell command used to invoke vampire-serial
//
// Returns true if all values pass, false otherwise.
//------------------------------------------------------------------------------
bool final_values_test(const std::string dir, const std::vector<expected_value_t> values, const std::string executable){

   // record where we started so we can return after the test
   std::string path = std::filesystem::current_path();

   // fixed-width label for aligned terminal output of first value
   std::stringstream test_name;
   test_name << "Testing " << values[0].label << " for " << dir;
   std::cout << std::setw(60) << std::left << test_name.str() << " : " << std::flush;

   // enter the test data directory
   if( !vt::chdir(path+"/data/"+dir) ) return false;

   // run vampire; non-zero exit code means something went wrong
   int vmp = vt::system(executable);
   if( vmp != 0){
      std::cerr << "Error running vampire. Returning as failed test." << std::endl;
      vt::remove({"output", "log"});
      vt::chdir(path);
      return false;
   }

   std::ifstream ifile;
   ifile.open("output");

   // scan through the entire output file, keeping only the last non-header
   // line. Header lines start with '#'; data lines start with a number.
   std::string line;
   std::string last_line;
   while(getline(ifile, line)){
      if(!line.empty() && line[0] != '#') last_line = line;
   }

   // clean up simulation output before returning (keeps the repo clean)
   vt::remove({"output", "log"});

   // return to the directory we started in
   if( !vt::chdir(path) ) return false;

   // parse all columns on the last data line as doubles
   std::stringstream ss(last_line);
   std::vector<double> columns;
   double column_value = 0.0;
   while(ss >> column_value) columns.push_back(column_value);

   // check each requested value
   bool pass = true;
   for(size_t i = 0; i < values.size(); ++i){

      // print fixed-width label for subsequent values (first already printed before running)
      if(i > 0){
         std::stringstream value_name;
         value_name << "Testing " << values[i].label << " for " << dir;
         std::cout << std::setw(60) << std::left << value_name.str() << " : ";
      }

      // check requested column exists in output
      if(values[i].column < 0 || values[i].column >= static_cast<int>(columns.size())){
         std::cout << "FAIL | column " << values[i].column << " not found in output" << std::endl;
         pass = false;
         continue;
      }

      // pass if |value - expected| <= abs_tol + tol*|expected|
      const double value = columns[values[i].column];
      const double expected = values[i].expected;
      if( std::abs(value - expected) <= values[i].abs_tol + values[i].tol * std::abs(expected) ){
         std::cout << "OK" << std::endl;
      }
      else{
         std::cout << "FAIL | expected: " << values[i].expected << "\tobtained: " << value << std::endl;
         pass = false;
      }

   }

   return pass;

}

//------------------------------------------------------------------------------
// final_value_test
//
// Runs a vampire simulation in the given subdirectory, reads the last scalar
// value written to the output file, and checks it matches 'expected' to
// within a relative tolerance 'tol' (i.e. |value/expected - 1| < tol).
//
// Parameters:
//   dir        - path relative to developers/test/integration/data/
//   expected   - the value the last output line should converge to
//   tol        - relative tolerance for the comparison (e.g. 1e-4)
//   executable - full shell command used to invoke vampire-serial
//
// Returns true if the test passes, false otherwise.
//------------------------------------------------------------------------------
bool final_value_test(const std::string dir, double expected, double tol, const std::string executable){

   return final_values_test(dir, { {"convergence", 0, expected, tol} }, executable);

}
