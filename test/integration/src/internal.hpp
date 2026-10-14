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
// include file only once
#pragma once

// C++ standard library headers
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// module headers
#include "internal.hpp"

namespace vt{
   int system(std::string const &cmd);
   bool chdir(const std::string path);
   void remove(const std::vector<std::string> files);
}

//------------------------------------------------------------------------------
// Test functions
//------------------------------------------------------------------------------
bool energy_test(std::string dir, double result, std::string executable, std::string energy_type);
bool integrator_test(const std::string dir, double time, double rx, double ry, double rz, const std::string executable);
bool material_atoms_test(const std::string dir, int n1, int n2, int n3, int n4, const std::string executable);
bool final_value_test(const std::string dir, double expected, double tol, const std::string executable);

// expected value of a single output column for final_values_test
struct expected_value_t{
   std::string label; // description printed in test output
   int column;        // 0-based column index in output file
   double expected;   // expected value
   double tol;        // relative tolerance
};
bool final_values_test(const std::string dir, const std::vector<expected_value_t> values, const std::string executable);
