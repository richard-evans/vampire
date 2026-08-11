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
#include <iostream>
#include <vector>

// Vampire headers
#include "create/internal.hpp"
#include "random.hpp"

// include header for test functions
#include "voronoi_bimodal_seeds_test.hpp"

namespace ut{
   namespace create{

//------------------------------------------------------------------------------
// No two placed seeds should ever be closer than the sum of the smallest
// radius each population's lognormal draw can produce after clamping (see
// generate_bimodal_voronoi_seeds(), r_min_clamp = 0.1*mean_diameter). This is
// a weaker, but always-valid, lower bound on the true per-seed exclusion
// distance, since individual drawn radii are not exposed by the function.
//------------------------------------------------------------------------------
int test_no_collisions(const bool verbose){

   int ec = 0;

   const double std_diameter = 50.0;   // Angstroms
   const double small_diameter = 20.0; // Angstroms
   const double variance = 0.15;

   const double min_r_std = 0.1*std_diameter;
   const double min_r_small = 0.1*small_diameter;

   const struct { double domain; double fraction; } cases[] = {
      { 200.0, 0.1 }, { 200.0, 0.3 }, { 200.0, 0.5 }, { 300.0, 0.9 }
   };

   for(size_t c=0; c<sizeof(cases)/sizeof(cases[0]); c++){

      mtrandom::grnd.seed(123456);

      std::vector<std::vector<double> > grain_coord_array;
      std::vector<bool> is_small_grain;

      ::create::internal::generate_bimodal_voronoi_seeds(
         grain_coord_array, is_small_grain,
         cases[c].domain, cases[c].domain,
         std_diameter, variance, small_diameter, variance, cases[c].fraction);

      const size_t n = grain_coord_array.size();
      for(size_t i=0; i<n; i++){
         const double min_ri = is_small_grain[i] ? min_r_small : min_r_std;
         for(size_t j=i+1; j<n; j++){
            const double min_rj = is_small_grain[j] ? min_r_small : min_r_std;
            const double dx = grain_coord_array[i][0] - grain_coord_array[j][0];
            const double dy = grain_coord_array[i][1] - grain_coord_array[j][1];
            const double dist = std::sqrt(dx*dx + dy*dy);
            const double min_allowed = min_ri + min_rj;
            if(dist < min_allowed - 1.0e-6){
               if(verbose) std::cout << "FAIL: seeds " << i << " and " << j << " are closer ("
                                      << dist << ") than the minimum possible exclusion distance ("
                                      << min_allowed << ") for fraction=" << cases[c].fraction << std::endl;
               ec++;
            }
         }
      }
   }

   return ec;
}

//------------------------------------------------------------------------------
// Two calls with the same mtrandom::voronoi_seed-derived RNG state must
// produce identical seed lists.
//------------------------------------------------------------------------------
int test_determinism(const bool verbose){

   int ec = 0;

   std::vector<std::vector<double> > coords_a, coords_b;
   std::vector<bool> small_a, small_b;

   mtrandom::grnd.seed(987654);
   ::create::internal::generate_bimodal_voronoi_seeds(
      coords_a, small_a, 250.0, 250.0, 50.0, 0.15, 20.0, 0.15, 0.3);

   mtrandom::grnd.seed(987654);
   ::create::internal::generate_bimodal_voronoi_seeds(
      coords_b, small_b, 250.0, 250.0, 50.0, 0.15, 20.0, 0.15, 0.3);

   if(coords_a.size() != coords_b.size()){
      if(verbose) std::cout << "FAIL: determinism - grain counts differ (" << coords_a.size()
                             << " vs " << coords_b.size() << ")" << std::endl;
      ec++;
      return ec;
   }

   for(size_t i=0; i<coords_a.size(); i++){
      if(coords_a[i][0] != coords_b[i][0] || coords_a[i][1] != coords_b[i][1] || small_a[i] != small_b[i]){
         if(verbose) std::cout << "FAIL: determinism - seed " << i << " differs between identically-seeded runs" << std::endl;
         ec++;
      }
   }

   return ec;
}

//------------------------------------------------------------------------------
// For a range of small_fraction values on a reasonably large domain, both
// populations should be non-empty and the realised fraction should be
// reasonably close to the request (loose statistical check, not exact - see
// manual entry for create:voronoi-bimodal-grains).
//------------------------------------------------------------------------------
int test_both_populations_present(const bool verbose){

   int ec = 0;

   const double fractions[] = { 0.1, 0.3, 0.5, 0.9 };
   const double tolerance = 0.15;

   for(size_t i=0; i<sizeof(fractions)/sizeof(fractions[0]); i++){

      mtrandom::grnd.seed(42);

      std::vector<std::vector<double> > grain_coord_array;
      std::vector<bool> is_small_grain;

      ::create::internal::generate_bimodal_voronoi_seeds(
         grain_coord_array, is_small_grain, 400.0, 400.0, 50.0, 0.15, 20.0, 0.15, fractions[i]);

      int placed_small = 0;
      for(size_t g=0; g<is_small_grain.size(); g++) if(is_small_grain[g]) placed_small++;
      const int placed_std = int(is_small_grain.size()) - placed_small;

      if(placed_small == 0){
         if(verbose) std::cout << "FAIL: no small grains placed for fraction=" << fractions[i] << std::endl;
         ec++;
      }
      if(placed_std == 0){
         if(verbose) std::cout << "FAIL: no standard grains placed for fraction=" << fractions[i] << std::endl;
         ec++;
      }

      if(!grain_coord_array.empty()){
         const double realised = double(placed_small)/double(grain_coord_array.size());
         if(std::fabs(realised - fractions[i]) > tolerance){
            if(verbose) std::cout << "FAIL: realised fraction " << realised << " too far from requested "
                                   << fractions[i] << " (tolerance " << tolerance << ")" << std::endl;
            ec++;
         }
      }
   }

   return ec;
}

//------------------------------------------------------------------------------
// Degenerate cases: small_fraction=0.0 places only standard-sized seeds, and
// a domain too small to fit even one grain terminates via the internal
// failure limit without hanging or crashing.
//------------------------------------------------------------------------------
int test_degenerate_cases(const bool verbose){

   int ec = 0;

   // small_fraction = 0.0 -> only standard seeds
   {
      mtrandom::grnd.seed(1);

      std::vector<std::vector<double> > grain_coord_array;
      std::vector<bool> is_small_grain;

      ::create::internal::generate_bimodal_voronoi_seeds(
         grain_coord_array, is_small_grain, 200.0, 200.0, 50.0, 0.15, 20.0, 0.15, 0.0);

      if(grain_coord_array.empty()){
         if(verbose) std::cout << "FAIL: small_fraction=0.0 placed no grains at all" << std::endl;
         ec++;
      }
      for(size_t g=0; g<is_small_grain.size(); g++){
         if(is_small_grain[g]){
            if(verbose) std::cout << "FAIL: small_fraction=0.0 placed a small grain at index " << g << std::endl;
            ec++;
            break;
         }
      }
   }

   // domain too small to fit even one grain comfortably still terminates
   {
      mtrandom::grnd.seed(2);

      std::vector<std::vector<double> > grain_coord_array;
      std::vector<bool> is_small_grain;

      ::create::internal::generate_bimodal_voronoi_seeds(
         grain_coord_array, is_small_grain, 0.0, 0.0, 50.0, 0.15, 20.0, 0.15, 0.3);

      // function returning at all (rather than hanging) is the test; a
      // degenerate domain should still place at least the first seed
      if(grain_coord_array.empty()){
         if(verbose) std::cout << "FAIL: degenerate zero-size domain placed no seeds at all" << std::endl;
         ec++;
      }
   }

   return ec;
}

//------------------------------------------------------------------------------
// Function to test create module functions
//------------------------------------------------------------------------------
int test_voronoi_bimodal_seeds(const bool verbose){

   if(verbose) std::cout << "Testing create::internal::generate_bimodal_voronoi_seeds()" << std::endl;

   int ec = 0;

   ec += test_no_collisions(verbose);
   ec += test_determinism(verbose);
   ec += test_both_populations_present(verbose);
   ec += test_degenerate_cases(verbose);

   return ec;

}

}
}
