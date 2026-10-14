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
// Usage: ./integration_tests [--extended | --all] [--mpi N]
//
//   (no option)  run standard tests only (fast, run on every change)
//   --extended   run extended tests only (long running statistical tests)
//   --all        run standard and extended tests
//   --mpi N      run vampire-parallel on N processes instead of vampire-serial
//------------------------------------------------------------------------------
int main(int argc, char* argv[]){

   // determine which test sets to run
   bool standard = true;
   bool extended = false;
   int mpi_processes = 0; // 0 = serial
   for(int arg = 1; arg < argc; ++arg){
      const std::string option = argv[arg];
      if(option == "--extended"){ standard = false; extended = true; }
      else if(option == "--all"){ standard = true;  extended = true; }
      else if(option == "--mpi" && arg + 1 < argc && std::atoi(argv[arg+1]) > 0){
         mpi_processes = std::atoi(argv[arg+1]);
         ++arg;
      }
      else{
         std::cerr << "Unknown option " << option << std::endl;
         std::cerr << "Usage: " << argv[0] << " [--extended | --all] [--mpi N]" << std::endl;
         return 1;
      }
   }

   std::cout << "---------------------------------------------------------------------" << std::endl;
   std::cout << "    Running system test suite for vampire code";
   if(mpi_processes > 0) std::cout << " (mpi, np = " << mpi_processes << ")";
   std::cout << std::endl;
   std::cout << "---------------------------------------------------------------------" << std::endl;

   // Get vampire path and command
   //std::string path = std::filesystem::current_path();
   std::filesystem::path wd=std::filesystem::current_path();
   std::string path_string = std::filesystem::current_path().parent_path().parent_path();

   std::string exe = path_string+"/vampire-serial 1>/dev/null";
   if(mpi_processes > 0){
      exe = "mpirun -np " + std::to_string(mpi_processes) + " " + path_string + "/vampire-parallel 1>/dev/null";
   }

   //std::cout << exe << std::endl;

   //return 0;

   //---------------------------------------------------------------------------
   // Run tests
   //---------------------------------------------------------------------------
   unsigned int fail = 0;

   //---------------------------------------------------------------------------
   // Standard tests
   //---------------------------------------------------------------------------
   if(standard){

      // Exchange energy tests
      if( !energy_test("crystals/sc" , -3.0e-17, exe, "exchange" ) ) fail += 1;
      if( !energy_test("crystals/fcc", -2.4e-16, exe, "exchange" ) ) fail += 1;
      if( !energy_test("crystals/bcc", -8.0e-18, exe, "exchange" ) ) fail += 1;

      // Anisotropy energy tests
      if( !energy_test("anisotropy/uniaxial", -1.0e-20, exe, "anisotropy" ) ) fail += 1;

      // DMI energy tests
      if( !energy_test("energy/dmi", -2.0e-18, exe, "dmi" ) ) fail += 1;

      // Integrator tests
      if( !integrator_test("dynamics/heun",    9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( !integrator_test("dynamics/midpoint",9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( !integrator_test("dynamics/rkmk2",   9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( !integrator_test("dynamics/rkmk4",   9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( !integrator_test("dynamics/dm",      9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;

      // Program tests
      // effective-damping: single spin rotated 120 degrees from z relaxing in 1T field along z
      if( !integrator_test("dynamics/effective-damping", 9.75e-11, -0.552854, 0.159192, 0.817931, exe ) ) fail += 1;

      // Structure tests
      if( !material_atoms_test("structure/core-shell", 3474, 485, 0, 0, exe ) ) fail += 1;
      if( !material_atoms_test("structure/multilayer", 1000, 1000, 0, 0, exe ) ) fail += 1;

      // Monte Carlo equilibration test
      if( !final_value_test("montecarlo/equilibrium", 1.0, 1.0e-4, exe ) ) fail += 1;

      // Constants tests
      // gamma_SI: single spin precession — spin trajectory at t=9.74e-11 s with alpha=0.001
      if( !integrator_test("constants/gamma", 9.74e-11, -0.127636, -0.991673, 0.0171491, exe ) ) fail += 1;
      // muB: Zeeman energy for 125 spins of 1 muB in 1T field = -125*muB = -1.15925125979e-21 J
      if( !energy_test("constants/muB", -1.15925125979e-21, exe, "muB" ) ) fail += 1;
      // kB: single spin MC with fixed seed; final spin direction is deterministic and depends on kB
      if( !integrator_test("constants/kB", 9.74e-11, 0.411154, -0.829297, 0.378442, exe ) ) fail += 1;

   }

   //---------------------------------------------------------------------------
   // Extended tests
   //
   // Long running tests (several seconds or more each), typically comparing
   // thermal averages against exact results, with tolerances set by the
   // statistical noise.
   //---------------------------------------------------------------------------
   if(extended){

      // Statistics tests
      // non-interacting spins (mu_s = 2 muB, B = 1 T along z, T = 1 K) with exact Langevin results
      // where x = mu_s B / kB T = 1.343428 and L(x) = coth(x) - 1/x:
      //    specific heat (kB per spin)        C/kB  = 1 - x^2/sinh^2(x)                 = 0.433945
      //    transverse susceptibility (1/T)    chi_x = chi_y = (mu_s/kB T) L(x)/x          = 0.401776
      //    longitudinal susceptibility (1/T)  chi_z = (mu_s/kB T) (1/x^2 - 1/sinh^2(x))  = 0.323013
      if( !final_values_test("statistics/paramagnet", { {"specific heat",    0, 0.433945, 5.0e-2},
                                                        {"susceptibility x", 1, 0.401776, 5.0e-2},
                                                        {"susceptibility y", 2, 0.401776, 5.0e-2},
                                                        {"susceptibility z", 3, 0.323013, 5.0e-2} }, exe ) ) fail += 1;

   }

   // Summary
   std::cout << "---------------------------------------------------------------------" << std::endl;
   if(fail >0){
      std::cout << "Failed " << fail << " tests : OVERALL FAIL" << std::endl;
   }
   else{
      std::cout << "Failed " << fail << " tests : OVERALL PASS" <<	std::endl;
   }
   std::cout << "---------------------------------------------------------------------" << std::endl;

   return fail > 0 ? 1 : 0;

}
