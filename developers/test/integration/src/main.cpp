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
   // the vampire root is the first directory above this one that contains hdr/
   std::filesystem::path root = std::filesystem::current_path();
   while(!std::filesystem::is_directory(root/"hdr") && root != root.parent_path()) root = root.parent_path();
   std::string path_string = root.string();

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
   unsigned int skipped = 0;

   // tests on a single spin cannot be decomposed over several mpi processes, so are skipped in mpi mode
   auto serial_only = [&](const std::string dir){
      if(mpi_processes == 0) return true;
      std::stringstream test_name;
      test_name << "Skipping " << dir;
      std::cout << std::setw(60) << std::left << test_name.str() << " : SKIPPED (single spin, serial only)" << std::endl;
      skipped += 1;
      return false;
   };

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

      // Anisotropy energy tests for each anisotropy term: five non-interacting layers of
      // 100 atoms with spins along [001], [100], [111], d4 = (0.3,-0.5,0.8) and d5 = (-0.6,0.2,0.4)
      // (normalised); output:material-total-energy gives the energy of each layer
      // uniaxial-second-order: E = -K (S.e)^2, e = [111], for 100 spins along each direction
      if( !final_values_test("anisotropy/uniaxial-second-order", { {"E[001]", 0, -3.333333333333334e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[100]", 1, -3.333333333333334e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[111]", 2, -1.000000000000000e-21, 1.0e-10, 1.0e-30},
                                                                   {"E(d4)", 3, -1.224489795918368e-22, 1.0e-10, 1.0e-30},
                                                                   {"E(d5)", 4, 0.0, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // uniaxial-fourth-order: E = -K [(S.e)^4 - 6/7 (S.e)^2], e = [111], for 100 spins along each direction
      if( !final_values_test("anisotropy/uniaxial-fourth-order", { {"E[001]", 0, 1.746031746031746e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[100]", 1, 1.746031746031746e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[111]", 2, -1.428571428571435e-22, 1.0e-10, 1.0e-30},
                                                                   {"E(d4)", 3, 8.996251561849234e-23, 1.0e-10, 1.0e-30},
                                                                   {"E(d5)", 4, 0.0, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // uniaxial-sixth-order: E = -K [(S.e)^6 - 15/11 (S.e)^4 + 5/11 (S.e)^2], e = [111], for 100 spins along each direction
      if( !final_values_test("anisotropy/uniaxial-sixth-order", { {"E[001]", 0, -3.703703703703701e-23, 1.0e-10, 1.0e-30},
                                                                  {"E[100]", 1, -3.703703703703701e-23, 1.0e-10, 1.0e-30},
                                                                  {"E[111]", 2, -9.090909090909143e-23, 1.0e-10, 1.0e-30},
                                                                  {"E(d4)", 3, -3.704857051676831e-23, 1.0e-10, 1.0e-30},
                                                                  {"E(d5)", 4, 0.0, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // cubic-fourth-order: E = -K/2 (Sx^4 + Sy^4 + Sz^4), for 100 spins along each direction
      if( !final_values_test("anisotropy/cubic-fourth-order", { {"E[001]", 0, -5.000000000000000e-22, 1.0e-10, 1.0e-30},
                                                                {"E[100]", 1, -5.000000000000000e-22, 1.0e-10, 1.0e-30},
                                                                {"E[111]", 2, -1.666666666666668e-22, 1.0e-10, 1.0e-30},
                                                                {"E(d4)", 3, -2.500000000000001e-22, 1.0e-10, 1.0e-30},
                                                                {"E(d5)", 4, -2.499999999999999e-22, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // cubic-sixth-order: E = K Sx^2 Sy^2 Sz^2, for 100 spins along each direction
      if( !final_values_test("anisotropy/cubic-sixth-order", { {"E[001]", 0, 0.0, 1.0e-10, 1.0e-30},
                                                               {"E[100]", 1, 0.0, 1.0e-10, 1.0e-30},
                                                               {"E[111]", 2, 3.703703703703706e-23, 1.0e-10, 1.0e-30},
                                                               {"E(d4)", 3, 1.529974755416535e-23, 1.0e-10, 1.0e-30},
                                                               {"E(d5)", 4, 1.311953352769679e-23, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // cubic-fourth-order-rotated: E = -K/2 sum_i (S.e_i)^4, e1 = [110], e2 = [-110], e3 = e1 x e2, for 100 spins along each direction
      if( !final_values_test("anisotropy/cubic-fourth-order-rotated", { {"E[001]", 0, -4.999999999999996e-22, 1.0e-10, 1.0e-30},
                                                                        {"E[100]", 1, -2.499999999999999e-22, 1.0e-10, 1.0e-30},
                                                                        {"E[111]", 2, -2.777777777777778e-22, 1.0e-10, 1.0e-30},
                                                                        {"E(d4)", 3, -2.667638483965012e-22, 1.0e-10, 1.0e-30},
                                                                        {"E(d5)", 4, -2.142857142857142e-22, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // biaxial-fourth-order: E = -K/2 [(S.u1)^4 + (S.u2)^4], u1 = [110], u2 = [1-10], for 100 spins along each direction
      if( !final_values_test("anisotropy/biaxial-fourth-order", { {"E[001]", 0, 0.0, 1.0e-10, 1.0e-30},
                                                                  {"E[100]", 1, -2.499999999999999e-22, 1.0e-10, 1.0e-30},
                                                                  {"E[111]", 2, -2.222222222222222e-22, 1.0e-10, 1.0e-30},
                                                                  {"E(d4)", 3, -5.351936693044559e-23, 1.0e-10, 1.0e-30},
                                                                  {"E(d5)", 4, -1.734693877551020e-22, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // triaxial-second-order: E = -(kx Sx^2 + ky Sy^2 + kz Sz^2), k = (1,2,-3) x 1e-23 J, for 100 spins along each direction
      if( !final_values_test("anisotropy/triaxial-second-order", { {"E[001]", 0, 3.000000000000000e-21, 1.0e-10, 1.0e-30},
                                                                   {"E[100]", 1, -9.999999999999999e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[111]", 2, 0.0, 1.0e-10, 1.0e-30},
                                                                   {"E(d4)", 3, 1.357142857142857e-21, 1.0e-10, 1.0e-30},
                                                                   {"E(d5)", 4, 7.142857142857171e-23, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // triaxial-fourth-order: E = -sum_i k_i (S_i^4 - 30/35 S_i^2), k = (1,2,-3) x 1e-23 J, for 100 spins along each direction
      if( !final_values_test("anisotropy/triaxial-fourth-order", { {"E[001]", 0, 4.285714285714287e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[100]", 1, -1.428571428571429e-22, 1.0e-10, 1.0e-30},
                                                                   {"E[111]", 2, -3.673419846319648e-38, 1.0e-10, 1.0e-30},
                                                                   {"E(d4)", 3, -2.238650562265698e-23, 1.0e-10, 1.0e-30},
                                                                   {"E(d5)", 4, -2.397959183673470e-22, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;

      // DMI energy tests
      if( !energy_test("energy/dmi", -2.0e-18, exe, "dmi" ) ) fail += 1;

      // Exchange energy tests on a chain of alternating A and B spins (chain.ucf) along general
      // directions S_A = (0.3,-0.5,0.8) and S_B = (0.5,0.4,0.7) (normalised)
      // isotropic exchange: E = -N J (S_A.S_B), J = 1e-21 J, for N = 2000 spins
      if( !final_values_test("exchange/isotropic", { {"exchange energy", 0, -1.086090160499898e-18, 1.0e-10} }, exe ) ) fail += 1;
      // vectorial exchange: E = -N sum_k J_kk S_A,k S_B,k, J = (1,2,3) x 1e-21 J, for N = 2000 spins
      if( !final_values_test("exchange/vectorial", { {"exchange energy", 0, -3.045311626499714e-18, 1.0e-10} }, exe ) ) fail += 1;
      // tensorial exchange: E = -N S_A^T J S_B for a general 3x3 J (with J^T for B->A bonds), for N = 2000 spins
      if( !final_values_test("exchange/tensorial", { {"exchange energy", 0, -2.606616385199755e-18, 1.0e-10} }, exe ) ) fail += 1;
      // biquadratic exchange: E = -N Jbq (S_A.S_B)^2, Jbq = 1e-21 J, bilinear J = 0, for N = 2000 spins
      if( !final_values_test("exchange/biquadratic", { {"exchange energy", 0, -5.897959183673468e-19, 1.0e-10} }, exe ) ) fail += 1;
      // kitaev exchange: E = -N J (S_A.S_B) + N K S_A,x S_B,x for bonds along x, J = 1e-21 J, K = 0.5e-21 J,
      // for N = 2000 spins (nearest neighbour bonds only, all within the Kitaev cutoff)
      if( !final_values_test("exchange/kitaev", { {"exchange energy", 0, -9.263710192499127e-19, 1.0e-10} }, exe ) ) fail += 1;

      // Integrator tests
      if( serial_only("dynamics/heun") && !integrator_test("dynamics/heun",    9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( serial_only("dynamics/midpoint") && !integrator_test("dynamics/midpoint",9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( serial_only("dynamics/rkmk2") && !integrator_test("dynamics/rkmk2",   9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( serial_only("dynamics/rkmk4") && !integrator_test("dynamics/rkmk4",   9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;
      if( serial_only("dynamics/dm") && !integrator_test("dynamics/dm",      9.74e-11, -0.103927, -0.338608, 0.935171, exe ) ) fail += 1;

      // Integrator convergence order: error at t = 20 ps for dt = 100 fs and 50 fs vs the exact damped precession
      if( serial_only("dynamics/convergence-llg-heun") && !convergence_test("dynamics/convergence-llg-heun", 0.1, 1.0, 2.0, 0.15, exe ) ) fail += 1;
      if( serial_only("dynamics/convergence-llg-midpoint") && !convergence_test("dynamics/convergence-llg-midpoint", 0.1, 1.0, 2.0, 0.15, exe ) ) fail += 1;
      if( serial_only("dynamics/convergence-llg-dm") && !convergence_test("dynamics/convergence-llg-dm", 0.1, 1.0, 2.0, 0.15, exe ) ) fail += 1;
      if( serial_only("dynamics/convergence-llg-rkmk2") && !convergence_test("dynamics/convergence-llg-rkmk2", 0.1, 1.0, 2.0, 0.15, exe ) ) fail += 1;
      if( serial_only("dynamics/convergence-llg-rkmk4") && !convergence_test("dynamics/convergence-llg-rkmk4", 0.1, 1.0, 4.0, 0.15, exe ) ) fail += 1;

      // Program tests
      // effective-damping: single spin rotated 120 degrees from z relaxing in 1T field along z
      if( serial_only("dynamics/effective-damping") && !integrator_test("dynamics/effective-damping", 9.75e-11, -0.552854, 0.159192, 0.817931, exe ) ) fail += 1;

      // Structure tests
      if( !material_atoms_test("structure/core-shell", 3474, 485, 0, 0, exe ) ) fail += 1;
      if( !material_atoms_test("structure/multilayer", 1000, 1000, 0, 0, exe ) ) fail += 1;

      // Monte Carlo equilibration test
      if( !final_value_test("montecarlo/equilibrium", 1.0, 1.0e-4, exe ) ) fail += 1;

      // Checkpoint tests: a run saved part way and continued must match an uninterrupted run exactly
      if( !checkpoint_test("checkpoint/llg-heun", exe ) ) fail += 1;
      if( !checkpoint_test("checkpoint/monte-carlo", exe ) ) fail += 1;

      // Constants tests
      // gamma_SI: single spin precession — spin trajectory at t=9.74e-11 s with alpha=0.001
      if( serial_only("constants/gamma") && !integrator_test("constants/gamma", 9.74e-11, -0.127636, -0.991673, 0.0171491, exe ) ) fail += 1;
      // muB: Zeeman energy for 125 spins of 1 muB in 1T field = -125*muB = -1.15925125979e-21 J
      if( !energy_test("constants/muB", -1.15925125979e-21, exe, "muB" ) ) fail += 1;

      // Zeeman energy: E = -mu_s B (S.h), mu_s = 2 muB, B = 1.5 T, h = (1,2,2)/3, for 100 spins along
      // [001], [100], [111], d4 = (0.3,-0.5,0.8) and d5 = (-0.6,0.2,0.4) (normalised)
      if( !final_values_test("energy/zeeman", { {"E[001]", 0, -1.854802015660000e-21, 1.0e-10, 1.0e-30},
                                                {"E[100]", 1, -9.274010078300000e-22, 1.0e-10, 1.0e-30},
                                                {"E[111]", 2, -2.677176107586904e-21, 1.0e-10, 1.0e-30},
                                                {"E(d4)", 3, -8.431348390917834e-22, 1.0e-10, 1.0e-30},
                                                {"E(d5)", 4, -7.435750353104395e-22, 1.0e-10, 1.0e-30} }, exe ) ) fail += 1;
      // kB: single spin MC with fixed seed; final spin direction is deterministic and depends on kB
      if( serial_only("constants/kB") && !integrator_test("constants/kB", 9.74e-11, 0.411154, -0.829297, 0.378442, exe ) ) fail += 1;

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
      std::cout << "Failed " << fail << " tests";
      if(skipped > 0) std::cout << ", skipped " << skipped;
      std::cout << " : OVERALL FAIL" << std::endl;
   }
   else{
      std::cout << "Failed " << fail << " tests";
      if(skipped > 0) std::cout << ", skipped " << skipped;
      std::cout << " : OVERALL PASS" << std::endl;
   }
   std::cout << "---------------------------------------------------------------------" << std::endl;

   return fail > 0 ? 1 : 0;

}
