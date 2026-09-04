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
#include "grains/internal.hpp"

// include header for test functions
#include "grain_shape_test.hpp"

namespace ut{
   namespace grains{

      typedef ::grains::internal::point2_t point2_t;
      typedef ::grains::internal::polygon_t polygon_t;

namespace{

   //--------------------------------------------------------------------
   // Local helpers (deliberately independent of grain_geometry_test.cpp's
   // own copies - each unit test file is self-contained, per this suite's
   // existing convention)
   //--------------------------------------------------------------------
   inline bool nearly(double a, double b, double tol = 1.0e-9){
      return std::fabs(a - b) < tol;
   }

   // wraps an angle difference to (-pi, pi]
   inline double angle_diff(double a, double b){
      double d = std::fmod(a - b + M_PI, 2.0*M_PI);
      if(d < 0.0) d += 2.0*M_PI;
      return d - M_PI;
   }

   polygon_t rectangle_poly(double xmin, double ymin, double xmax, double ymax){
      return ::grains::internal::rectangle(xmin, ymin, xmax, ymax);
   }

   // an irregular (but convex) pentagon, centred roughly on the origin,
   // used wherever a genuinely non-regular cell is needed (containment/
   // area-bound cases, where the point is that faceting must clip it)
   polygon_t irregular_pentagon(){
      polygon_t poly;
      poly.push_back(point2_t( 5.0, -3.0));
      poly.push_back(point2_t( 6.0,  2.0));
      poly.push_back(point2_t( 1.0,  6.0));
      poly.push_back(point2_t(-4.0,  3.0));
      poly.push_back(point2_t(-2.0, -4.0));
      return poly;
   }

   // true if p lies inside (or on the boundary of, within tol) the convex
   // CCW polygon poly
   bool point_in_convex(const polygon_t& poly, point2_t p, double tol = 1.0e-7){
      const size_t n = poly.size();
      if(n < 3) return false;
      for(size_t i=0; i<n; i++){
         const point2_t& a = poly[i];
         const point2_t& b = poly[(i+1)%n];
         const double ex = b.x - a.x, ey = b.y - a.y;
         const double cross = ex*(p.y - a.y) - ey*(p.x - a.x);
         if(cross < -tol) return false;
      }
      return true;
   }

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Case a: strength = 0 (or num_facets <= 0) leaves the cell bit-identical.
//------------------------------------------------------------------------------
int test_strength_zero_unchanged(const bool verbose){

   int ec = 0;

   const polygon_t poly = irregular_pentagon();
   const point2_t centre(0.3, -0.2); // deliberately off-centre from the pentagon's own centroid

   const polygon_t r1 = ::grains::internal::apply_wulff_facets(poly, centre, 6, 0.37, 1.3, 0.0, 100.0);
   const polygon_t r2 = ::grains::internal::apply_wulff_facets(poly, centre, 0, 0.37, 1.3, 1.0, 100.0);

   if(r1.size() != poly.size()){
      if(verbose) std::cout << "FAIL: strength=0 changed vertex count (" << r1.size() << " vs " << poly.size() << ")" << std::endl;
      ec++;
   }
   else{
      for(size_t v=0; v<poly.size(); v++){
         if(r1[v].x != poly[v].x || r1[v].y != poly[v].y){
            if(verbose) std::cout << "FAIL: strength=0 vertex " << v << " changed: (" << r1[v].x << "," << r1[v].y
                                   << ") vs original (" << poly[v].x << "," << poly[v].y << ")" << std::endl;
            ec++;
         }
      }
   }

   if(r2.size() != poly.size()){
      if(verbose) std::cout << "FAIL: num_facets=0 changed vertex count (" << r2.size() << " vs " << poly.size() << ")" << std::endl;
      ec++;
   }
   else{
      for(size_t v=0; v<poly.size(); v++){
         if(r2[v].x != poly[v].x || r2[v].y != poly[v].y){
            if(verbose) std::cout << "FAIL: num_facets=0 vertex " << v << " changed" << std::endl;
            ec++;
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case b: strength = 1, m = 4, isotropic gamma, applied to a large cell -> a
// square of the requested (target_area) area, rotated by theta0. The "large
// cell" (a big axis-aligned square, far larger than the target) exists so
// the pure Wulff square fits entirely inside it and the clip removes
// nothing - see internal.hpp's own note on why target_area is an
// explicit parameter, independent of the cell's own area, rather than
// derived from it (equal-area containment is only possible when the cell
// already IS the target shape).
//------------------------------------------------------------------------------
int test_isotropic_square(const bool verbose){

   int ec = 0;

   const point2_t centre(0.0, 0.0);
   const polygon_t large_cell = rectangle_poly(-100.0, -100.0, 100.0, 100.0);

   const double target_area = 400.0; // side 20, apothem s = sqrt(400/4) = 10
   const double theta0 = 0.4;
   const double s = std::sqrt(target_area/4.0);

   const polygon_t result = ::grains::internal::apply_wulff_facets(large_cell, centre, 4, theta0, 1.0, 1.0, target_area);

   if(result.size() != 4){
      if(verbose) std::cout << "FAIL: isotropic square facet count = " << result.size() << ", expected 4" << std::endl;
      return ec+1;
   }

   const double area = ::grains::internal::polygon_area(result);
   if(!nearly(area, target_area, 1.0e-6)){
      if(verbose) std::cout << "FAIL: isotropic square area = " << area << ", expected " << target_area << std::endl;
      ec++;
   }

   // vertices of a square facetted at distance s from centre in 4 directions
   // theta0+90k lie at angle theta0+45+90k, radius s*sqrt(2) (see grain_shape_test.cpp
   // derivation: intersection of two facet planes 90 degrees apart, both at
   // distance s, sits at the angular bisector, radius s/cos(45deg))
   const double expected_radius = s/std::cos(M_PI/4.0);
   for(int k=0; k<4; k++){
      const double expected_angle = theta0 + M_PI/4.0 + k*M_PI/2.0;
      const point2_t expected(centre.x + expected_radius*std::cos(expected_angle),
                               centre.y + expected_radius*std::sin(expected_angle));

      bool found = false;
      for(size_t v=0; v<result.size(); v++){
         if(nearly(result[v].x, expected.x, 1.0e-6) && nearly(result[v].y, expected.y, 1.0e-6)){ found = true; break; }
      }
      if(!found){
         if(verbose) std::cout << "FAIL: isotropic square missing expected vertex ("
                                << expected.x << "," << expected.y << ")" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case c: m = 6 -> a regular hexagon; vertex angles correct to 1e-9.
//------------------------------------------------------------------------------
int test_regular_hexagon(const bool verbose){

   int ec = 0;

   const point2_t centre(1.5, -0.5); // off-origin, to pin that centre (not the origin) is what matters
   const polygon_t large_cell = rectangle_poly(centre.x-100.0, centre.y-100.0, centre.x+100.0, centre.y+100.0);

   const double target_area = 3.0*std::sqrt(3.0)*10.0*10.0/2.0; // regular hexagon, apothem s=10: A = 6*s^2*tan(30deg) = 2*sqrt(3)*s^2... see below
   // For a regular m-gon of apothem s: A = m*s^2*tan(pi/m). For m=6: A = 6*s^2*tan(30deg) = 2*sqrt(3)*s^2.
   const double s = 10.0;
   const double expected_area = 6.0*s*s*std::tan(M_PI/6.0);
   const double theta0 = -0.9;

   const polygon_t result = ::grains::internal::apply_wulff_facets(large_cell, centre, 6, theta0, 1.0, 1.0, expected_area);

   if(result.size() != 6){
      if(verbose) std::cout << "FAIL: hexagon facet count = " << result.size() << ", expected 6" << std::endl;
      return ec+1;
   }

   const double area = ::grains::internal::polygon_area(result);
   if(!nearly(area, expected_area, 1.0e-6)){
      if(verbose) std::cout << "FAIL: hexagon area = " << area << ", expected " << expected_area << std::endl;
      ec++;
   }

   const double expected_radius = s/std::cos(M_PI/6.0);
   for(int k=0; k<6; k++){
      const double expected_angle = theta0 + M_PI/6.0 + k*M_PI/3.0;

      bool found = false;
      for(size_t v=0; v<result.size(); v++){
         const double dx = result[v].x - centre.x;
         const double dy = result[v].y - centre.y;
         const double radius = std::sqrt(dx*dx + dy*dy);
         const double angle = std::atan2(dy, dx);
         if(nearly(radius, expected_radius, 1.0e-6) && std::fabs(angle_diff(angle, expected_angle)) < 1.0e-9){
            found = true;
            break;
         }
      }
      if(!found){
         if(verbose) std::cout << "FAIL: hexagon missing expected vertex at angle " << expected_angle
                                << ", radius " << expected_radius << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case d: the faceted polygon is always contained in the original cell -
// swept over facet counts, anisotropy and strength, including
// configurations where the Wulff target does NOT fit (an irregular,
// deliberately non-roomy cell), which is exactly where containment is a
// real constraint rather than a triviality.
//------------------------------------------------------------------------------
int test_containment(const bool verbose){

   int ec = 0;

   const polygon_t poly = irregular_pentagon();
   const point2_t centre(0.0, 0.0); // pentagon roughly centred on the origin
   const double cell_area = ::grains::internal::polygon_area(poly);

   const int facet_counts[] = {4, 6, 8};
   const double strengths[] = {0.0, 0.3, 0.6, 1.0};
   const double anisotropies[] = {1.0, 0.4, 2.5};
   const double target_areas[] = {0.3*cell_area, cell_area, 3.0*cell_area}; // undersized, matched, oversized targets

   for(int mi=0; mi<3; mi++){
      for(int si=0; si<4; si++){
         for(int ai=0; ai<3; ai++){
            for(int ti=0; ti<3; ti++){

               const polygon_t result = ::grains::internal::apply_wulff_facets(
                  poly, centre, facet_counts[mi], 0.55, anisotropies[ai], strengths[si], target_areas[ti]);

               for(size_t v=0; v<result.size(); v++){
                  if(!point_in_convex(poly, result[v])){
                     if(verbose) std::cout << "FAIL: faceted vertex (" << result[v].x << "," << result[v].y
                                            << ") lies outside the original cell (m=" << facet_counts[mi]
                                            << ", strength=" << strengths[si] << ", anisotropy=" << anisotropies[ai]
                                            << ", target_area=" << target_areas[ti] << ")" << std::endl;
                     ec++;
                  }
               }

            }
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case e: faceted area <= cell area for all lambda (same sweep as case d).
//------------------------------------------------------------------------------
int test_area_bound(const bool verbose){

   int ec = 0;

   const polygon_t poly = irregular_pentagon();
   const point2_t centre(0.0, 0.0);
   const double cell_area = ::grains::internal::polygon_area(poly);

   const int facet_counts[] = {4, 6, 8};
   const double strengths[] = {0.0, 0.3, 0.6, 1.0};
   const double anisotropies[] = {1.0, 0.4, 2.5};
   const double target_areas[] = {0.3*cell_area, cell_area, 3.0*cell_area};

   for(int mi=0; mi<3; mi++){
      for(int si=0; si<4; si++){
         for(int ai=0; ai<3; ai++){
            for(int ti=0; ti<3; ti++){

               const polygon_t result = ::grains::internal::apply_wulff_facets(
                  poly, centre, facet_counts[mi], 0.55, anisotropies[ai], strengths[si], target_areas[ti]);

               const double area = ::grains::internal::polygon_area(result);
               if(area > cell_area + 1.0e-6){
                  if(verbose) std::cout << "FAIL: faceted area " << area << " exceeds cell area " << cell_area
                                        << " (m=" << facet_counts[mi] << ", strength=" << strengths[si]
                                        << ", anisotropy=" << anisotropies[ai] << ", target_area=" << target_areas[ti]
                                        << ")" << std::endl;
                  ec++;
               }

            }
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// generate_grain_facet_orientations() pins each of its three modes
// directly: random orientation stays within [0,2*pi/m) (facets have
// m-fold rotational symmetry, so this is the full range of distinct
// orientations), fixed is exact for every grain, and textured scatters
// around the fixed angle rather than being exactly fixed.
//------------------------------------------------------------------------------
int test_generate_orientations(const bool verbose){

   int ec = 0;

   MTRand rng;
   rng.seed(20260812);

   // random: every value lies in [0, 2*pi/m)
   {
      const int m = 6;
      const std::vector<double> theta = ::grains::internal::generate_grain_facet_orientations(
         2000, m, ::grains::internal::facet_orientation_random, 0.0, 0.0, rng);
      const double period = 2.0*M_PI/double(m);
      for(size_t i=0; i<theta.size(); i++){
         if(theta[i] < 0.0 || theta[i] >= period){
            if(verbose) std::cout << "FAIL: random orientation " << theta[i] << " outside [0," << period << ")" << std::endl;
            ec++;
            break;
         }
      }
   }

   // fixed: every grain gets exactly the requested angle
   {
      const double angle_degrees = 37.0;
      const std::vector<double> theta = ::grains::internal::generate_grain_facet_orientations(
         50, 4, ::grains::internal::facet_orientation_fixed, angle_degrees, 0.0, rng);
      const double expected = angle_degrees*M_PI/180.0;
      for(size_t i=0; i<theta.size(); i++){
         if(theta[i] != expected){
            if(verbose) std::cout << "FAIL: fixed orientation " << theta[i] << " != " << expected << std::endl;
            ec++;
            break;
         }
      }
   }

   // textured: scatters around the fixed angle (not all identical, unlike "fixed")
   {
      const double angle_degrees = 10.0;
      const std::vector<double> theta = ::grains::internal::generate_grain_facet_orientations(
         2000, 4, ::grains::internal::facet_orientation_textured, angle_degrees, 15.0, rng);
      const double expected = angle_degrees*M_PI/180.0;
      bool any_different = false;
      double sum = 0.0;
      for(size_t i=0; i<theta.size(); i++){
         sum += theta[i];
         if(!nearly(theta[i], expected, 1.0e-9)) any_different = true;
      }
      if(!any_different){
         if(verbose) std::cout << "FAIL: textured orientation did not scatter around the fixed angle" << std::endl;
         ec++;
      }
      const double mean = sum/double(theta.size());
      if(!nearly(mean, expected, 0.05)){ // 0.05 rad ~ 3deg, generous for 2000 draws at 15deg sd
         if(verbose) std::cout << "FAIL: textured orientation mean " << mean << " far from fixed angle " << expected << std::endl;
         ec++;
      }
   }

   // num_facets<=0 (faceting off): well-defined all-zero array
   {
      const std::vector<double> theta = ::grains::internal::generate_grain_facet_orientations(
         10, 0, ::grains::internal::facet_orientation_random, 0.0, 0.0, rng);
      if(theta.size() != 10){
         if(verbose) std::cout << "FAIL: num_facets<=0 orientation array size = " << theta.size() << ", expected 10" << std::endl;
         ec++;
      }
      for(size_t i=0; i<theta.size(); i++){
         if(theta[i] != 0.0){
            if(verbose) std::cout << "FAIL: num_facets<=0 orientation[" << i << "] = " << theta[i] << ", expected 0.0" << std::endl;
            ec++;
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Function to test create module Wulff faceting
//------------------------------------------------------------------------------
int test_grain_shape(const bool verbose){

   if(verbose) std::cout << "Testing grains::internal:: Wulff faceting" << std::endl;

   int ec = 0;

   ec += test_strength_zero_unchanged(verbose);
   ec += test_isotropic_square(verbose);
   ec += test_regular_hexagon(verbose);
   ec += test_containment(verbose);
   ec += test_area_bound(verbose);
   ec += test_generate_orientations(verbose);

   return ec;

}

   }
}
