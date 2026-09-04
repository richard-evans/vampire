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
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire headers
#include "grains/internal.hpp"

// include header for test functions
#include "grain_packing_test.hpp"

namespace ut{
   namespace grains{

      typedef ::grains::internal::point2_t point2_t;

namespace{

   // n discs of diameter d, sized so their target areas sum to roughly the
   // domain area (the same balance choose_grain_count() aims for in
   // production) - a representative, non-dilute packing problem.
   std::vector<double> uniform_diameters(int n, double d){
      return std::vector<double>(n, d);
   }

} // end of anonymous namespace

//------------------------------------------------------------------------------
// A representative near-full-density packing (20 discs of diameter 10 in a
// ~40x40 domain, area fraction ~1) settles to a bounded overlap rather than
// diverging or leaving discs wildly interpenetrating. Circles cannot tile a
// plane at 100% area fraction, so some residual overlap at this density is
// expected and left for fit_grain_weights() to correct - this only pins
// that relaxation stays within a sane bound rather than requiring exact
// convergence.
//------------------------------------------------------------------------------
int test_grain_packing_bounded_overlap(const bool verbose){

   int ec = 0;

   const int n = 20;
   const double d = 10.0;
   const double domain = std::sqrt(n * M_PI/4.0 * d*d);

   const std::vector<double> diameters = uniform_diameters(n, d);

   MTRand rng;
   rng.seed(20260826);

   const ::grains::internal::grain_packing_result_t result =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, false, rng);

   if(result.sites.size() != size_t(n)){
      if(verbose) std::cout << "FAIL: grain packing returned " << result.sites.size()
                             << " sites, expected " << n << std::endl;
      return ec+1;
   }

   if(!std::isfinite(result.achieved_max_overlap) || result.achieved_max_overlap < 0.0 || result.achieved_max_overlap > 0.5){
      if(verbose) std::cout << "FAIL: grain packing achieved max overlap fraction " << result.achieved_max_overlap
                             << " after " << result.iterations << " iterations, expected in [0, 0.5]" << std::endl;
      ec++;
   }

   for(size_t i = 0; i < result.sites.size(); i++){
      if(!std::isfinite(result.sites[i].x) || !std::isfinite(result.sites[i].y)){
         if(verbose) std::cout << "FAIL: site " << i << " has a non-finite position" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Reproducibility: the same seed gives an identical packed site list.
//------------------------------------------------------------------------------
int test_grain_packing_reproducibility(const bool verbose){

   int ec = 0;

   const std::vector<double> diameters = uniform_diameters(15, 8.0);
   const double domain = 35.0;

   MTRand rng;

   rng.seed(13371337);
   const ::grains::internal::grain_packing_result_t a =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, false, rng);

   rng.seed(13371337);
   const ::grains::internal::grain_packing_result_t b =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, false, rng);

   if(a.sites.size() != b.sites.size()){
      if(verbose) std::cout << "FAIL: packing reproducibility - site counts differ (" << a.sites.size()
                             << " vs " << b.sites.size() << ")" << std::endl;
      return ec+1;
   }

   for(size_t i = 0; i < a.sites.size(); i++){
      if(a.sites[i].x != b.sites[i].x || a.sites[i].y != b.sites[i].y){
         if(verbose) std::cout << "FAIL: packing reproducibility - site " << i << " differs between identically-seeded runs" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Non-periodic mode: wall repulsion must keep every disc centre within the
// requested [0,domain] bounds (never truncated past the box, never placed
// outside it).
//------------------------------------------------------------------------------
int test_grain_packing_nonperiodic_bounds(const bool verbose){

   int ec = 0;

   const std::vector<double> diameters = uniform_diameters(25, 6.0);
   const double domain = 25.0;

   MTRand rng;
   rng.seed(90210);

   const ::grains::internal::grain_packing_result_t result =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, false, rng);

   for(size_t i = 0; i < result.sites.size(); i++){
      const point2_t& p = result.sites[i];
      if(p.x < 0.0 || p.x > domain || p.y < 0.0 || p.y > domain){
         if(verbose) std::cout << "FAIL: non-periodic site " << i << " at (" << p.x << "," << p.y
                                << ") lies outside [0," << domain << "]" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Periodic mode applies no wall force, so - unlike non-periodic mode -
// discs are free to relax with centres arbitrarily close to the domain
// edge (they wrap, rather than being pushed inward). Comparing the closest
// any site gets to an edge under periodic vs. non-periodic relaxation for
// the same seed/diameters/domain pins that periodic mode carries no
// wall-driven push-away bias: non-periodic's minimum edge distance should
// stay close to a full disc radius (the wall force's approximate
// stand-off), periodic's should not.
//------------------------------------------------------------------------------
int test_grain_packing_periodic_no_wall_bias(const bool verbose){

   int ec = 0;

   // Kept at a modest area fraction (~0.3, not the near-full-density
   // packing test_grain_packing_bounded_overlap uses): at high density the
   // wall's soft repulsion can legitimately be pushed through by pressure
   // from the interior, since it is a spring rather than a hard barrier
   // (per design - "pushed inward ... never truncated or reflected"), which
   // would blur the periodic/non-periodic distinction this test is pinning.
   const std::vector<double> diameters = uniform_diameters(30, 5.0);
   const double domain = 45.0;
   const double mean_radius = 2.5;

   auto min_edge_distance = [&](const ::grains::internal::grain_packing_result_t& r){
      double worst = domain;
      for(size_t i = 0; i < r.sites.size(); i++){
         const point2_t& p = r.sites[i];
         worst = std::min(worst, std::min(p.x, domain - p.x));
         worst = std::min(worst, std::min(p.y, domain - p.y));
      }
      return worst;
   };

   MTRand rng;

   rng.seed(555555);
   const ::grains::internal::grain_packing_result_t periodic =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, true, rng);

   rng.seed(555555);
   const ::grains::internal::grain_packing_result_t nonperiodic =
      ::grains::internal::pack_grain_seeds(diameters, domain, domain, false, rng);

   const double d_periodic = min_edge_distance(periodic);
   const double d_nonperiodic = min_edge_distance(nonperiodic);

   // non-periodic's wall force keeps every centre a non-trivial distance
   // clear of the edge; periodic has no such force, so its closest
   // approach should be comfortably smaller. The bound is well below a
   // full radius (rather than ~mean_radius) because the wall is a soft
   // spring, not a hard barrier - under packing pressure from interior
   // neighbours it can still be pushed some way through, just consistently
   // less far than an edge with no wall force at all.
   if(d_nonperiodic < 0.2*mean_radius){
      if(verbose) std::cout << "FAIL: non-periodic minimum edge distance " << d_nonperiodic
                             << " is smaller than expected for an active wall force (radius " << mean_radius << ")" << std::endl;
      ec++;
   }

   if(d_periodic >= d_nonperiodic){
      if(verbose) std::cout << "FAIL: periodic minimum edge distance " << d_periodic
                             << " is not smaller than non-periodic's " << d_nonperiodic
                             << " - periodic packing should show no wall-driven push-away bias" << std::endl;
      ec++;
   }

   return ec;

}

//------------------------------------------------------------------------------
int test_grain_packing(const bool verbose){

   if(verbose) std::cout << "Testing grains::internal:: physical disc-packing seed placement" << std::endl;

   int ec = 0;

   ec += test_grain_packing_bounded_overlap(verbose);
   ec += test_grain_packing_reproducibility(verbose);
   ec += test_grain_packing_nonperiodic_bounds(verbose);
   ec += test_grain_packing_periodic_no_wall_bias(verbose);

   return ec;

}

   }
}
