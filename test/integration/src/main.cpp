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
   if( !exchange_test("crystals/sc" , -3.0e-17, exe ) ) fail += 1;
   if( !exchange_test("crystals/fcc", -2.4e-16, exe ) ) fail += 1;

   // Integrator tests
   if( !integrator_test("dynamics/heun",-0.106813,-0.337996,0.935067, exe ) ) fail += 1;

   // Structure tests
   if( !material_atoms_test("structure/core-shell", 3474, 485, 0, 0, exe ) ) fail += 1;

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
