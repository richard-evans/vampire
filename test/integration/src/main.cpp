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

int main(){

   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "    Running system test suite for vampire code" << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;

   // Get vampire path and command
   //std::string path = std::filesystem::current_path();
   std::filesystem::path wd=std::filesystem::current_path();
   std::string path_string = std::filesystem::current_path().parent_path().parent_path();

   std::string exe = path_string+"/vampire-serial 1>/dev/null";

   //std::cout << exe << std::endl;

   //return 0;

   //---------------------------------------------------------------------------
   // Run tests
   //---------------------------------------------------------------------------
   unsigned int fail = 0;

   // Exchange energy tests
   if( !energy_test("crystals/sc" , -3.0e-17, exe, "exchange" ) ) fail += 1;
   if( !energy_test("crystals/fcc", -2.4e-16, exe, "exchange" ) ) fail += 1;
   if( !energy_test("crystals/bcc", -8.0e-18, exe, "exchange" ) ) fail += 1;

   // Anisotropy energy tests
   if( !energy_test("anisotropy/uniaxial", -1.0e-20, exe, "anisotropy" ) ) fail += 1;

   // Integrator tests
   if( !integrator_test("dynamics/heun",    -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
   if( !integrator_test("dynamics/midpoint",-0.0849683,-0.339983, 0.936585, exe ) ) fail += 1;

   // Structure tests
   if( !material_atoms_test("structure/core-shell", 3474, 485, 0, 0, exe ) ) fail += 1;
   if( !material_atoms_test("structure/multilayer", 1000, 1000, 0, 0, exe ) ) fail += 1;

   // Monte Carlo equilibration test
   if( !final_value_test("montecarlo/equilibrium", 1.0, 1.0e-4, exe ) ) fail += 1;

   // Constants tests
   // gamma_SI: single spin precession — spin trajectory at t=9.74e-11 s with alpha=0.001
   if( !integrator_test("constants/gamma", -0.127636, -0.991673, 0.0171491, exe ) ) fail += 1;
   // muB: Zeeman energy for 125 spins of 1 muB in 1T field = -125*muB = -1.15925125979e-21 J
   if( !energy_test("constants/muB", -1.15925125979e-21, exe, "muB" ) ) fail += 1;
   // kB: single spin MC with fixed seed; final spin direction is deterministic and depends on kB
   if( !integrator_test("constants/kB", 0.411154, -0.829297, 0.378442, exe ) ) fail += 1;

   // Summary
   std::cout << "---------------------------------------------------------------------" << std::endl;
   if(fail >0){
      std::cout << "Failed " << fail << " tests : OVERALL FAIL" << std::endl;
   }
   else{
      std::cout << "Failed " << fail << " tests : OVERALL PASS" <<	std::endl;
   }
   std::cout << "---------------------------------------------------------------------" << std::endl;

   return 0;

}
