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

// module headers
#include "internal.hpp"

namespace vt{
   int system(std::string const &cmd);
   bool chdir(const std::string path);
}

//------------------------------------------------------------------------------
// Test functions
//------------------------------------------------------------------------------
bool exchange_test(std::string dir, double result, std::string executable);
bool integrator_test(const std::string dir, double rx, double ry, double rz, const std::string executable);
bool material_atoms_test(const std::string dir, int n1, int n2, int n3, int n4, const std::string executable);
