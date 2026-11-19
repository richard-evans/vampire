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
#include <iostream>
#include <list>
#include <vector>

// Vampire headers
#include "create.hpp"
#include "material.hpp"

// include header for test functions
#include "create_tests.hpp"

// Material globals read by create::sorted_core_shell_materials(); in vampire they are defined in
// src/main/initialise_variables.cpp, which the unit tests do not link
namespace mp{
   int num_materials = 1;
   std::vector<materials_t> material(1);
}

namespace ut{

//------------------------------------------------------------------------------
// Core-shell material order: decreasing radius, equal radii in decreasing material
// index, and a strict weak ordering (std::list::sort is undefined otherwise)
//------------------------------------------------------------------------------
int core_shell_order_tests(const bool verbose){

   int error_count = 0;

   // irreflexive: an element never comes before itself
   for(double r : {0.0, 0.5, 1.0}){
      const create::core_radius_t a = {2, r};
      if(create::compare_radius(a, a)){
         std::cout << "   compare_radius(a,a) is true for radius " << r << std::endl;
         error_count++;
      }
   }

   // asymmetric for equal radii: the higher material index comes first
   const create::core_radius_t m0 = {0, 1.0};
   const create::core_radius_t m1 = {1, 1.0};
   if(!create::compare_radius(m1, m0) || create::compare_radius(m0, m1)){
      std::cout << "   equal radii are not ordered by decreasing material index" << std::endl;
      error_count++;
   }

   // sorted list: radii 0.5, 1, 1, 0.5, 1 for materials 0-4 give 4 2 1 3 0 (the order of the
   // earlier comparison function with libstdc++ and libc++, kept so material assignment is unchanged)
   const std::vector<double> radii = {0.5, 1.0, 1.0, 0.5, 1.0};
   mp::num_materials = radii.size();
   mp::material.resize(radii.size());
   for(size_t mat = 0; mat < radii.size(); mat++) mp::material[mat].core_shell_size = radii[mat];
   const std::vector<int> expected = {4, 2, 1, 3, 0};
   std::vector<int> obtained;
   for(const create::core_radius_t& c : create::sorted_core_shell_materials()) obtained.push_back(c.mat);
   if(obtained != expected){
      std::cout << "   sorted core-shell materials:";
      for(int mat : obtained) std::cout << " " << mat;
      std::cout << " (expected 4 2 1 3 0)" << std::endl;
      error_count++;
   }

   // all radii equal (the default core-shell-size): decreasing material index
   for(size_t mat = 0; mat < radii.size(); mat++) mp::material[mat].core_shell_size = 1.0;
   obtained.clear();
   for(const create::core_radius_t& c : create::sorted_core_shell_materials()) obtained.push_back(c.mat);
   if(obtained != std::vector<int>({4, 3, 2, 1, 0})){
      std::cout << "   materials with equal core-shell sizes are not in decreasing index order" << std::endl;
      error_count++;
   }

   mp::num_materials = 1;
   mp::material.resize(1);

   if(verbose) std::cout << "   core-shell order     : " << (error_count == 0 ? "PASS" : "FAIL") << std::endl;
   return error_count;

}

//------------------------------------------------------------------------------
// Function to test create module functions
//------------------------------------------------------------------------------
int create_tests(const bool verbose){

   if(verbose) std::cout << "Testing create module" << std::endl;

   int error_count = 0;

   error_count += core_shell_order_tests(verbose);

   if(verbose) std::cout <<          "================================" << std::endl;
   if(error_count == 0) std::cout << " create               : PASS " << std::endl;
   else std::cout <<                 " create               : FAIL " << error_count << std::endl;
   if(verbose) std::cout <<          "================================" << std::endl;

   return error_count;

}

}
