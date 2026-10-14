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

// Vampire headers
#include "grains/internal.hpp"
#include "random.hpp"

// include header for test functions
#include "grain_geometry_test.hpp"

namespace ut{
   namespace grains{

      typedef ::grains::internal::point2_t point2_t;
      typedef ::grains::internal::polygon_t polygon_t;

namespace{

   //--------------------------------------------------------------------
   // Local helpers
   //--------------------------------------------------------------------
   inline bool nearly(double a, double b, double tol = 1.0e-9){
      return std::fabs(a - b) < tol;
   }

   polygon_t unit_square(){
      return ::grains::internal::rectangle(0.0, 0.0, 1.0, 1.0);
   }

   // regular polygon of m sides, circumradius r, centred at (cx,cy), CCW
   polygon_t regular_polygon(int m, double r, double cx = 0.0, double cy = 0.0){
      polygon_t poly;
      for(int k=0; k<m; k++){
         const double theta = 2.0*M_PI*double(k)/double(m);
         poly.push_back(point2_t(cx + r*std::cos(theta), cy + r*std::sin(theta)));
      }
      return poly;
   }

   // True if poly is star-shaped about centre with CCW winding: walking its
   // vertices in order, the polar angle about centre must be monotonically
   // non-decreasing and sweep out exactly one full turn. This is exactly
   // the property round_polygon()'s organic-rounding overload relies on
   // (a strictly positive, single-valued radius as a function of angle is
   // star-shaped and simple by construction), and is a strictly weaker
   // check than convexity - a star-shaped polygon can have inward-bulging
   // (locally concave) edges, which the roughness overload is expected to
   // produce, while still being guaranteed simple (non-self-intersecting).
   bool is_star_shaped_ccw(const polygon_t& poly, point2_t centre, double tol = 1.0e-9){
      const size_t n = poly.size();
      if(n < 3) return true;

      double prev_theta = std::atan2(poly[0].y-centre.y, poly[0].x-centre.x);
      double total_turn = 0.0;
      for(size_t i=1; i<=n; i++){
         const point2_t& p = poly[i%n];
         const double theta = std::atan2(p.y-centre.y, p.x-centre.x);
         double dtheta = theta - prev_theta;
         while(dtheta <= -M_PI) dtheta += 2.0*M_PI;
         while(dtheta >  M_PI) dtheta -= 2.0*M_PI;
         if(dtheta < -tol) return false;
         total_turn += dtheta;
         prev_theta = theta;
      }
      return std::fabs(total_turn - 2.0*M_PI) < 1.0e-6;
   }

   // cross-product sign test: true if poly is convex and CCW (or trivially
   // small, n<3)
   bool is_convex_ccw(const polygon_t& poly, double tol = 1.0e-9){
      const size_t n = poly.size();
      if(n < 3) return true;
      for(size_t i=0; i<n; i++){
         const point2_t& a = poly[i];
         const point2_t& b = poly[(i+1)%n];
         const point2_t& c = poly[(i+2)%n];
         const double ex1 = b.x - a.x, ey1 = b.y - a.y;
         const double ex2 = c.x - b.x, ey2 = c.y - b.y;
         const double cross = ex1*ey2 - ey1*ex2;
         if(cross < -tol) return false;
      }
      return true;
   }

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Cases a, b: polygon_area()
//------------------------------------------------------------------------------
int test_polygon_area(const bool verbose){

   int ec = 0;

   // a. unit square area = 1
   {
      const double area = ::grains::internal::polygon_area(unit_square());
      if(!nearly(area, 1.0)){
         if(verbose) std::cout << "FAIL: unit square area = " << area << ", expected 1.0" << std::endl;
         ec++;
      }
   }

   // a. regular hexagon of circumradius r: area = 3*sqrt(3)/2 * r^2
   {
      const double r = 2.0;
      const double expected = 3.0*std::sqrt(3.0)/2.0*r*r;
      const double area = ::grains::internal::polygon_area(regular_polygon(6, r));
      if(!nearly(area, expected, 1.0e-9)){
         if(verbose) std::cout << "FAIL: hexagon area = " << area << ", expected " << expected << std::endl;
         ec++;
      }
   }

   // b. signed area positive for CCW, negative for CW (reversed) input
   {
      polygon_t ccw = unit_square();
      polygon_t cw(ccw.rbegin(), ccw.rend());

      const double area_ccw = ::grains::internal::polygon_area(ccw);
      const double area_cw  = ::grains::internal::polygon_area(cw);

      if(area_ccw <= 0.0){
         if(verbose) std::cout << "FAIL: CCW square area not positive (" << area_ccw << ")" << std::endl;
         ec++;
      }
      if(area_cw >= 0.0){
         if(verbose) std::cout << "FAIL: CW square area not negative (" << area_cw << ")" << std::endl;
         ec++;
      }
      if(!nearly(area_ccw, -area_cw)){
         if(verbose) std::cout << "FAIL: CCW/CW areas not equal and opposite (" << area_ccw << " vs " << area_cw << ")" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Cases c, d, e, f, g: clip_halfplane()
//------------------------------------------------------------------------------
int test_clip_halfplane(const bool verbose){

   int ec = 0;

   // c. clip unit square by x <= 0.5 -> area 0.5, 4 vertices
   {
      const polygon_t clipped = ::grains::internal::clip_halfplane(unit_square(), 1.0, 0.0, 0.5);
      const double area = ::grains::internal::polygon_area(clipped);
      if(!nearly(area, 0.5)){
         if(verbose) std::cout << "FAIL: clip by x<=0.5 area = " << area << ", expected 0.5" << std::endl;
         ec++;
      }
      if(clipped.size() != 4){
         if(verbose) std::cout << "FAIL: clip by x<=0.5 vertex count = " << clipped.size() << ", expected 4" << std::endl;
         ec++;
      }
   }

   // d. clip unit square by the diagonal x+y <= 1 -> triangle, area 0.5, 3 vertices
   {
      const double inv_sqrt2 = 1.0/std::sqrt(2.0);
      const polygon_t clipped = ::grains::internal::clip_halfplane(unit_square(), inv_sqrt2, inv_sqrt2, inv_sqrt2);
      const double area = ::grains::internal::polygon_area(clipped);
      if(!nearly(area, 0.5, 1.0e-8)){
         if(verbose) std::cout << "FAIL: clip by diagonal area = " << area << ", expected 0.5" << std::endl;
         ec++;
      }
      if(clipped.size() != 3){
         if(verbose) std::cout << "FAIL: clip by diagonal vertex count = " << clipped.size() << ", expected 3" << std::endl;
         ec++;
      }
   }

   // e. clip entirely outside -> empty polygon, area 0, no crash
   {
      const polygon_t clipped = ::grains::internal::clip_halfplane(unit_square(), 1.0, 0.0, -1.0);
      if(!clipped.empty()){
         if(verbose) std::cout << "FAIL: clip entirely outside left " << clipped.size() << " vertices, expected empty" << std::endl;
         ec++;
      }
      if(!nearly(::grains::internal::polygon_area(clipped), 0.0)){
         if(verbose) std::cout << "FAIL: clip entirely outside area not zero" << std::endl;
         ec++;
      }
   }

   // f. clip entirely inside -> unchanged polygon
   {
      const polygon_t original = unit_square();
      const polygon_t clipped = ::grains::internal::clip_halfplane(original, 1.0, 0.0, 10.0);
      if(clipped.size() != original.size()){
         if(verbose) std::cout << "FAIL: clip entirely inside changed vertex count (" << clipped.size() << " vs " << original.size() << ")" << std::endl;
         ec++;
      }
      else{
         for(size_t i=0; i<original.size(); i++){
            if(!nearly(clipped[i].x, original[i].x) || !nearly(clipped[i].y, original[i].y)){
               if(verbose) std::cout << "FAIL: clip entirely inside moved vertex " << i << std::endl;
               ec++;
            }
         }
      }
   }

   // g. repeated clipping by the same plane is idempotent
   {
      const polygon_t once = ::grains::internal::clip_halfplane(unit_square(), 1.0, 0.3, 0.6);
      const polygon_t twice = ::grains::internal::clip_halfplane(once, 1.0, 0.3, 0.6);

      if(once.size() != twice.size()){
         if(verbose) std::cout << "FAIL: repeated clip not idempotent - vertex counts " << once.size() << " vs " << twice.size() << std::endl;
         ec++;
      }
      else{
         for(size_t i=0; i<once.size(); i++){
            if(!nearly(once[i].x, twice[i].x) || !nearly(once[i].y, twice[i].y)){
               if(verbose) std::cout << "FAIL: repeated clip not idempotent at vertex " << i << std::endl;
               ec++;
            }
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// ensure_ccw(): normalises vertex winding to counter-clockwise. Several
// downstream primitives (offset_inward, round_polygon, the facet clipper)
// compute an inward-facing normal directly from edge direction, which is
// only correct for CCW polygons; feeding them a clockwise polygon silently
// clips against outward-facing half-planes instead and collapses the shape.
//------------------------------------------------------------------------------
int test_ensure_ccw(const bool verbose){

   int ec = 0;

   const polygon_t ccw = unit_square();
   const polygon_t cw(ccw.rbegin(), ccw.rend());

   const polygon_t fixed = ::grains::internal::ensure_ccw(cw);
   if(fixed.size() != cw.size() || ::grains::internal::polygon_area(fixed) <= 0.0){
      if(verbose) std::cout << "FAIL: ensure_ccw did not fix CW winding (area = "
                             << ::grains::internal::polygon_area(fixed) << ")" << std::endl;
      ec++;
   }

   const polygon_t unchanged = ::grains::internal::ensure_ccw(ccw);
   if(unchanged.size() != ccw.size() || !nearly(::grains::internal::polygon_area(unchanged), ::grains::internal::polygon_area(ccw))){
      if(verbose) std::cout << "FAIL: ensure_ccw altered an already-CCW polygon" << std::endl;
      ec++;
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case h: offset_inward()
//------------------------------------------------------------------------------
int test_offset_inward(const bool verbose){

   int ec = 0;

   // h. offset_inward(square side s, delta) -> square of side s-2*delta
   {
      const double s = 4.0;
      const double delta = 0.75;
      const polygon_t square = ::grains::internal::rectangle(0.0, 0.0, s, s);
      const polygon_t offset = ::grains::internal::offset_inward(square, delta);

      const double expected_area = (s - 2.0*delta)*(s - 2.0*delta);
      const double area = ::grains::internal::polygon_area(offset);

      if(!nearly(area, expected_area, 1.0e-8)){
         if(verbose) std::cout << "FAIL: offset_inward square area = " << area << ", expected " << expected_area << std::endl;
         ec++;
      }
      if(offset.size() != 4){
         if(verbose) std::cout << "FAIL: offset_inward square vertex count = " << offset.size() << ", expected 4" << std::endl;
         ec++;
      }
   }

   // h. over-offsetting collapses to empty, not a self-intersecting polygon
   {
      const double s = 1.0;
      const polygon_t square = ::grains::internal::rectangle(0.0, 0.0, s, s);
      const polygon_t offset = ::grains::internal::offset_inward(square, 0.6); // > s/2

      if(!offset.empty()){
         if(verbose) std::cout << "FAIL: over-offset square left " << offset.size()
                                << " vertices (area " << ::grains::internal::polygon_area(offset)
                                << "), expected empty" << std::endl;
         ec++;
      }
   }

   // offset_inward on a regular hexagon reduces the apothem (the
   // centre-to-edge-midpoint distance) by exactly delta: moving every edge
   // inward along its own normal by delta is precisely what shrinks the
   // inradius/apothem by delta for a convex polygon, independent of shape.
   {
      const double r = 5.0;
      const double delta = 1.2;
      const polygon_t hex = regular_polygon(6, r);
      const double apothem0 = r*std::cos(M_PI/6.0);
      const double expected_apothem = apothem0 - delta;

      const polygon_t offset = ::grains::internal::offset_inward(hex, delta);

      if(offset.size() != 6){
         if(verbose) std::cout << "FAIL: offset_inward hexagon vertex count = " << offset.size() << ", expected 6" << std::endl;
         ec++;
      }
      else{
         for(size_t i=0; i<offset.size(); i++){
            const point2_t& p0 = offset[i];
            const point2_t& p1 = offset[(i+1)%offset.size()];
            const double mx = 0.5*(p0.x+p1.x);
            const double my = 0.5*(p0.y+p1.y);
            const double apothem = std::sqrt(mx*mx+my*my);
            if(!nearly(apothem, expected_apothem, 1.0e-8)){
               if(verbose) std::cout << "FAIL: offset_inward hexagon edge " << i << " apothem = " << apothem
                                      << ", expected " << expected_apothem << std::endl;
               ec++;
            }
         }
      }
   }

   // The gap between two offset NEIGHBOURING cells of different sizes,
   // sharing an edge, equals 2*delta regardless of which cell is bigger.
   // This is the defining property of offsetting each cell inward by a
   // fixed distance rather than by a size-dependent scale factor: a grain
   // boundary of constant physical width, independent of grain size.
   {
      const double delta = 0.8;

      // a large and a small rectangle sharing the edge x=0
      const polygon_t big   = ::grains::internal::rectangle(-10.0, -6.0, 0.0, 6.0);
      const polygon_t small = ::grains::internal::rectangle(0.0, -1.5, 2.0, 1.5);

      const polygon_t big_offset   = ::grains::internal::offset_inward(big, delta);
      const polygon_t small_offset = ::grains::internal::offset_inward(small, delta);

      double max_x_big = big_offset.empty() ? 0.0 : big_offset[0].x;
      for(size_t i=1; i<big_offset.size(); i++){
         if(big_offset[i].x > max_x_big) max_x_big = big_offset[i].x;
      }

      double min_x_small = small_offset.empty() ? 0.0 : small_offset[0].x;
      for(size_t i=1; i<small_offset.size(); i++){
         if(small_offset[i].x < min_x_small) min_x_small = small_offset[i].x;
      }

      const double gap = min_x_small - max_x_big;

      if(!nearly(gap, 2.0*delta, 1.0e-8)){
         if(verbose) std::cout << "FAIL: offset neighbour gap = " << gap << ", expected " << 2.0*delta
                                << " (independent of cell size)" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// round_polygon(): clips a polygon at a bisected radius so the retained
// area matches a requested fraction of the original area exactly.
//------------------------------------------------------------------------------
int test_round_polygon(const bool verbose){

   int ec = 0;

   // c. rounding a square to fraction f gives area f*A to 1e-6
   {
      const polygon_t square = ::grains::internal::rectangle(-2.0, -2.0, 2.0, 2.0);
      const double A = ::grains::internal::polygon_area(square);
      const double f = 0.6;

      const polygon_t rounded = ::grains::internal::round_polygon(square, point2_t(0.0, 0.0), f);
      const double area = ::grains::internal::polygon_area(rounded);

      if(!nearly(area, f*A, 1.0e-6)){
         if(verbose) std::cout << "FAIL: round_polygon square f=" << f << " area = " << area
                                << ", expected " << f*A << std::endl;
         ec++;
      }
   }

   // d. rounding with fraction 1.0 leaves the polygon unchanged
   {
      const polygon_t square = ::grains::internal::rectangle(-2.0, -2.0, 2.0, 2.0);
      const polygon_t rounded = ::grains::internal::round_polygon(square, point2_t(0.0, 0.0), 1.0);

      if(rounded.size() != square.size()){
         if(verbose) std::cout << "FAIL: round_polygon f=1.0 changed vertex count (" << rounded.size()
                                << " vs " << square.size() << ")" << std::endl;
         ec++;
      }
      else{
         for(size_t i=0; i<square.size(); i++){
            if(!nearly(rounded[i].x, square[i].x) || !nearly(rounded[i].y, square[i].y)){
               if(verbose) std::cout << "FAIL: round_polygon f=1.0 moved vertex " << i << std::endl;
               ec++;
            }
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// round_polygon()'s organic-rounding overload: rounds via a random per-
// facet radius instead of a single common one.
//------------------------------------------------------------------------------
int test_round_polygon_organic(const bool verbose){

   int ec = 0;

   // roughness_amplitude<=0.0 delegates to the plain overload, bit-identically
   {
      const polygon_t square = ::grains::internal::rectangle(-2.0, -2.0, 2.0, 2.0);
      const double f = 0.6;

      mtrandom::grnd.seed(20260826);
      const polygon_t plain = ::grains::internal::round_polygon(square, point2_t(0.0, 0.0), f);
      const polygon_t organic_off = ::grains::internal::round_polygon(square, point2_t(0.0, 0.0), f,
                                                                        0.0, 4, mtrandom::grnd);

      if(plain.size() != organic_off.size()){
         if(verbose) std::cout << "FAIL: round_polygon roughness=0 changed vertex count ("
                                << organic_off.size() << " vs " << plain.size() << ")" << std::endl;
         ec++;
      }
      else{
         for(size_t i=0; i<plain.size(); i++){
            if(!nearly(plain[i].x, organic_off[i].x) || !nearly(plain[i].y, organic_off[i].y)){
               if(verbose) std::cout << "FAIL: round_polygon roughness=0 is not bit-identical to the plain overload at vertex " << i << std::endl;
               ec++;
            }
         }
      }
   }

   // a positive roughness amplitude still hits the requested area fraction exactly
   {
      const polygon_t hex = regular_polygon(8, 5.0);
      const double A = ::grains::internal::polygon_area(hex);
      const double f = 0.5;

      mtrandom::grnd.seed(1);
      const polygon_t rounded = ::grains::internal::round_polygon(hex, point2_t(0.0, 0.0), f, 0.5, 4, mtrandom::grnd);
      const double area = ::grains::internal::polygon_area(rounded);

      if(!nearly(area, f*A, 1.0e-6)){
         if(verbose) std::cout << "FAIL: round_polygon organic f=" << f << " area = " << area
                                << ", expected " << f*A << std::endl;
         ec++;
      }
   }

   // the result stays a simple, star-shaped polygon (required by
   // vmath::point_in_polygon_scaled()'s radial-scaling assumption) across a
   // handful of random draws, at the maximum amplitude the interface parser
   // permits (create:grain-boundary-roughness is clamped to 0.6)
   {
      const polygon_t octagon = regular_polygon(8, 5.0);

      for(uint32_t seed=1; seed<=20; seed++){
         mtrandom::grnd.seed(seed);
         const polygon_t rounded = ::grains::internal::round_polygon(octagon, point2_t(0.0, 0.0), 0.7, 0.6, 5, mtrandom::grnd);

         if(!is_star_shaped_ccw(rounded, point2_t(0.0, 0.0))){
            if(verbose) std::cout << "FAIL: round_polygon organic result is not star-shaped/simple for seed " << seed << std::endl;
            ec++;
         }
         if(rounded.size() >= 3 && ::grains::internal::polygon_area(rounded) < -1.0e-9){
            if(verbose) std::cout << "FAIL: round_polygon organic result has negative area for seed " << seed << std::endl;
            ec++;
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case i: polygon_centroid()
//------------------------------------------------------------------------------
int test_polygon_centroid(const bool verbose){

   int ec = 0;

   // i. centroid of a symmetric polygon (square, off-origin) is its centre
   {
      const polygon_t square = ::grains::internal::rectangle(2.0, 3.0, 6.0, 9.0);
      const point2_t centroid = ::grains::internal::polygon_centroid(square);
      if(!nearly(centroid.x, 4.0) || !nearly(centroid.y, 6.0)){
         if(verbose) std::cout << "FAIL: square centroid = (" << centroid.x << "," << centroid.y
                                << "), expected (4,6)" << std::endl;
         ec++;
      }
   }

   // i. centroid of a regular polygon (off-origin) is its centre
   {
      const double cx = -3.0, cy = 5.0;
      const polygon_t hex = regular_polygon(6, 2.5, cx, cy);
      const point2_t centroid = ::grains::internal::polygon_centroid(hex);
      if(!nearly(centroid.x, cx, 1.0e-8) || !nearly(centroid.y, cy, 1.0e-8)){
         if(verbose) std::cout << "FAIL: hexagon centroid = (" << centroid.x << "," << centroid.y
                                << "), expected (" << cx << "," << cy << ")" << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Case j: clipping a convex polygon always yields a convex polygon
//------------------------------------------------------------------------------
int test_convexity_preserved(const bool verbose){

   int ec = 0;

   // an irregular (but convex) pentagon
   polygon_t pentagon;
   pentagon.push_back(point2_t(0.0, 0.0));
   pentagon.push_back(point2_t(4.0, -1.0));
   pentagon.push_back(point2_t(6.0, 2.0));
   pentagon.push_back(point2_t(3.0, 5.0));
   pentagon.push_back(point2_t(-1.0, 3.0));

   const polygon_t shapes[] = { unit_square(), regular_polygon(6, 3.0), pentagon };

   // a variety of clip planes, some through the interior, some missing it,
   // some through vertices/edges exactly
   const struct { double nx, ny, c; } planes[] = {
      { 1.0, 0.0,  0.5 },
      { 0.0, 1.0,  0.5 },
      { 1.0, 1.0,  1.0 },
      { 1.0, -1.0, 0.0 },
      { -1.0, 0.3, 0.2 },
      { 0.7, 0.7,  2.0 }
   };

   for(size_t s=0; s<sizeof(shapes)/sizeof(shapes[0]); s++){
      polygon_t poly = shapes[s];
      for(size_t p=0; p<sizeof(planes)/sizeof(planes[0]); p++){

         poly = ::grains::internal::clip_halfplane(poly, planes[p].nx, planes[p].ny, planes[p].c);

         if(!is_convex_ccw(poly)){
            if(verbose) std::cout << "FAIL: shape " << s << " lost convexity/CCW ordering after clip " << p << std::endl;
            ec++;
         }

         if(poly.size() >= 3 && ::grains::internal::polygon_area(poly) < -1.0e-9){
            if(verbose) std::cout << "FAIL: shape " << s << " has negative area after clip " << p << std::endl;
            ec++;
         }

         if(poly.empty()) break; // nothing left to keep clipping
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Function to test create module grain geometry primitives
//------------------------------------------------------------------------------
int test_grain_geometry(const bool verbose){

   if(verbose) std::cout << "Testing grains::internal grain geometry primitives" << std::endl;

   int ec = 0;

   ec += test_polygon_area(verbose);
   ec += test_clip_halfplane(verbose);
   ec += test_ensure_ccw(verbose);
   ec += test_offset_inward(verbose);
   ec += test_polygon_centroid(verbose);
   ec += test_convexity_preserved(verbose);
   ec += test_round_polygon(verbose);
   ec += test_round_polygon_organic(verbose);

   return ec;

}

}
}
