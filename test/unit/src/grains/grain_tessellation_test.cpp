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

// C++ standard library headers
#include <cmath>
#include <iostream>

// Vampire headers
#include "grains/internal.hpp"
#include "random.hpp"

// include header for test functions
#include "grain_tessellation_test.hpp"

namespace ut{
   namespace grains{

      typedef ::grains::internal::point2_t point2_t;
      typedef ::grains::internal::polygon_t polygon_t;

namespace{

   inline bool nearly(double a, double b, double tol = 1.0e-9){
      return std::fabs(a - b) < tol;
   }

   // cross-product-sign test: true if pt lies inside (or on the boundary
   // of) the CCW convex polygon poly
   bool point_in_convex_polygon(const polygon_t& poly, const point2_t& pt, double tol = 1.0e-9){
      const size_t n = poly.size();
      if(n < 3) return false;
      for(size_t i=0; i<n; i++){
         const point2_t& a = poly[i];
         const point2_t& b = poly[(i+1)%n];
         const double cross = (b.x-a.x)*(pt.y-a.y) - (b.y-a.y)*(pt.x-a.x);
         if(cross < -tol) return false;
      }
      return true;
   }

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Case a: 3x3 sites on a regular square lattice, uniform weights -> the
// interior cell is a square of the lattice spacing, to 1e-12. The diagonal
// neighbours sit exactly on the corners of that square (distance s/sqrt(2)
// from centre = half their bisector distance s*sqrt(2)/2), so they must not
// distort it.
//------------------------------------------------------------------------------
int test_square_lattice(const bool verbose){

   int ec = 0;

   const double s = 10.0; // lattice spacing

   std::vector<point2_t> sites;
   for(int gx=0; gx<3; gx++)
      for(int gy=0; gy<3; gy++)
         sites.push_back(point2_t(gx*s, gy*s));

   const std::vector<double> weights(sites.size(), 0.0); // uniform -> plain Voronoi

   const std::vector<polygon_t> cells = ::grains::internal::build_power_cells(
      sites, weights, -s, -s, 3.0*s, 3.0*s, false);

   // centre site is index 4 (gx=1,gy=1 -> position (s,s))
   const polygon_t& centre_cell = cells[4];

   const double area = ::grains::internal::polygon_area(centre_cell);
   if(!nearly(area, s*s, 1.0e-9)){
      if(verbose) std::cout << "FAIL: interior lattice cell area = " << area << ", expected " << s*s << std::endl;
      ec++;
   }

   for(size_t v=0; v<centre_cell.size(); v++){
      const double x = centre_cell[v].x - s; // relative to centre site
      const double y = centre_cell[v].y - s;
      if(!nearly(std::fabs(x), s*0.5, 1.0e-9) || !nearly(std::fabs(y), s*0.5, 1.0e-9)){
         if(verbose) std::cout << "FAIL: interior lattice cell vertex (" << x << "," << y
                                << ") is not a corner of the expected " << s << "x" << s << " square" << std::endl;
         ec++;
      }
   }

   if(centre_cell.size() != 4){
      if(verbose) std::cout << "FAIL: interior lattice cell has " << centre_cell.size()
                             << " vertices, expected 4" << std::endl;
      ec++;
   }

   return ec;

}

//------------------------------------------------------------------------------
// Two sites, weighted as touching discs (r1=60, r2=20, centre separation
// d=80=r1+r2), must produce a power-diagram boundary that passes exactly
// through the discs' point of tangency, i.e. at distance r1 from site 1.
// This is the defining property of the Laguerre/power diagram: unlike a
// plain Voronoi diagram, whose boundary is the perpendicular bisector of
// the two sites regardless of their size, the power diagram's boundary is
// the radical axis of the two weighted points (w_i = r_i^2), which is
// exactly the common tangent line for two touching circles. It is what
// lets each grain's cell contain exactly its own disc, whatever the size
// ratio between neighbours.
//------------------------------------------------------------------------------
int test_touching_discs(const bool verbose){

   int ec = 0;

   const double r1 = 60.0, r2 = 20.0, d = 80.0; // r1+r2 == d: touching

   const std::vector<point2_t> sites = { point2_t(0.0, 0.0), point2_t(d, 0.0) };
   const std::vector<double> weights = { r1*r1, r2*r2 };

   const std::vector<polygon_t> cells = ::grains::internal::build_power_cells(
      sites, weights, -100.0, -100.0, 180.0, 100.0, false);

   double max_x = -1.0e30;
   for(size_t v=0; v<cells[0].size(); v++) max_x = std::max(max_x, cells[0][v].x);

   if(!nearly(max_x, r1, 1.0e-9)){
      if(verbose) std::cout << "FAIL: touching-discs divider at x=" << max_x
                             << ", expected exactly r1=" << r1 << std::endl;
      ec++;
   }

   for(size_t v=0; v<cells[0].size(); v++){
      if(cells[0][v].x > r1 + 1.0e-9){
         if(verbose) std::cout << "FAIL: site-1 cell has a vertex at x=" << cells[0][v].x
                                << " beyond the tangent point r1=" << r1 << std::endl;
         ec++;
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Cases c, d, e, g on a shared random configuration: random sites (rejection-
// sampled to a minimum spacing) and random weights within a box.
//
// One property checked below is that every site lies inside its own cell.
// For a genuine power diagram this is NOT guaranteed for arbitrary weights:
// site i is excluded from its own cell ("orphaned") whenever a neighbour j
// satisfies |p_i-p_j|^2 < w_j - w_i, i.e. j's weight advantage outweighs
// the squared distance between them. Orphaning is correct Laguerre-diagram
// behaviour, not a defect, but it is a different property from the one this
// test is meant to isolate, so radii here are kept well under half the
// enforced minimum spacing, which guarantees |p_i-p_j|^2 >= min_spacing^2 >
// max_radius^2 >= w_j - w_i for every pair and rules orphaning out.
//------------------------------------------------------------------------------
int test_random_configuration(const bool verbose){

   int ec = 0;

   const double box = 100.0;
   const int n = 30;
   const double min_spacing = 8.0;
   const double max_radius = 2.0; // << min_spacing/2 - see note above

   mtrandom::grnd.seed(20260811);

   std::vector<point2_t> sites;
   std::vector<double> weights;

   int attempts = 0;
   while(int(sites.size()) < n && attempts < 100000){
      attempts++;
      const point2_t candidate(mtrandom::grnd()*box, mtrandom::grnd()*box);
      bool ok = true;
      for(size_t i=0; i<sites.size(); i++){
         const double dx = candidate.x - sites[i].x;
         const double dy = candidate.y - sites[i].y;
         if(dx*dx + dy*dy < min_spacing*min_spacing){ ok = false; break; }
      }
      if(ok){
         sites.push_back(candidate);
         const double r = mtrandom::grnd()*max_radius;
         weights.push_back(r*r);
      }
   }

   if(int(sites.size()) != n){
      if(verbose) std::cout << "FAIL: could not place " << n << " well-spaced random sites (placed "
                             << sites.size() << ")" << std::endl;
      ec++;
      return ec;
   }

   const std::vector<polygon_t> cells = ::grains::internal::build_power_cells(
      sites, weights, 0.0, 0.0, box, box, false);

   // c. partition of unity: sum of cell areas equals the box area
   {
      double total_area = 0.0;
      for(size_t i=0; i<cells.size(); i++) total_area += ::grains::internal::polygon_area(cells[i]);

      const double expected = box*box;
      if(!nearly(total_area, expected, 1.0e-6)){
         if(verbose) std::cout << "FAIL: partition of unity - total cell area = " << total_area
                                << ", expected " << expected << std::endl;
         ec++;
      }
   }

   for(size_t i=0; i<cells.size(); i++){

      const polygon_t& cell = cells[i];

      // e. no cell has zero or negative signed area, and no degenerate
      // 1-2 vertex artefact
      if(cell.size() == 1 || cell.size() == 2){
         if(verbose) std::cout << "FAIL: site " << i << " has a degenerate " << cell.size()
                                << "-vertex cell" << std::endl;
         ec++;
      }
      else if(cell.size() >= 3){
         const double area = ::grains::internal::polygon_area(cell);
         if(area < 1.0e-9){
            if(verbose) std::cout << "FAIL: site " << i << " has non-positive cell area " << area << std::endl;
            ec++;
         }

         // d. each site lies inside its own (non-empty) cell
         if(!point_in_convex_polygon(cell, sites[i])){
            if(verbose) std::cout << "FAIL: site " << i << " does not lie inside its own cell" << std::endl;
            ec++;
         }
      }

      // g. cells are clipped to the box exactly - no vertex outside
      for(size_t v=0; v<cell.size(); v++){
         if(cell[v].x < -1.0e-9 || cell[v].x > box+1.0e-9 || cell[v].y < -1.0e-9 || cell[v].y > box+1.0e-9){
            if(verbose) std::cout << "FAIL: site " << i << " has a vertex (" << cell[v].x << "," << cell[v].y
                                   << ") outside the domain [0," << box << "]^2" << std::endl;
            ec++;
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Under periodic boundaries the tessellated cells must still exactly
// partition the simulation box: their areas must sum to the box area, with
// no gaps or overlaps introduced by wrapping neighbours across the
// boundary. Uses the same rejection-sampled random configuration as
// test_random_configuration() above, with weights kept small relative to
// spacing so that no site is orphaned (see the note there).
//------------------------------------------------------------------------------
int test_periodic_partition_of_unity(const bool verbose){

   int ec = 0;

   const double box = 100.0;
   const int n = 20;
   const double min_spacing = 10.0;
   const double max_radius = 2.0; // << min_spacing/2

   mtrandom::grnd.seed(20260812);

   std::vector<point2_t> sites;
   std::vector<double> weights;

   int attempts = 0;
   while(int(sites.size()) < n && attempts < 100000){
      attempts++;
      const point2_t candidate(mtrandom::grnd()*box, mtrandom::grnd()*box);
      bool ok = true;
      for(size_t i=0; i<sites.size(); i++){
         const double dx = candidate.x - sites[i].x;
         const double dy = candidate.y - sites[i].y;
         if(dx*dx + dy*dy < min_spacing*min_spacing){ ok = false; break; }
      }
      if(ok){
         sites.push_back(candidate);
         const double r = mtrandom::grnd()*max_radius;
         weights.push_back(r*r);
      }
   }

   if(int(sites.size()) != n){
      if(verbose) std::cout << "FAIL: could not place " << n << " well-spaced random sites (placed "
                             << sites.size() << ")" << std::endl;
      ec++;
      return ec;
   }

   const std::vector<polygon_t> cells = ::grains::internal::build_power_cells(
      sites, weights, 0.0, 0.0, box, box, true);

   double total_area = 0.0;
   for(size_t i=0; i<cells.size(); i++) total_area += ::grains::internal::polygon_area(cells[i]);

   const double expected = box*box;
   if(!nearly(total_area, expected, 1.0e-6)){
      if(verbose) std::cout << "FAIL: periodic partition of unity - total cell area = " << total_area
                             << ", expected " << expected << std::endl;
      ec++;
   }

   return ec;

}

//------------------------------------------------------------------------------
// A periodic tessellation must be translation-covariant: shifting every
// site (and the domain itself) by one full box length is physically the
// same periodic configuration viewed in a differently-anchored frame, so
// every cell must come back translated by exactly the same amount. The
// whole coordinate frame is shifted together, rather than moving a single
// site outside [xmin,xmax] with the box held fixed, so that every site
// stays validly bounded within the domain passed to build_power_cells().
//------------------------------------------------------------------------------
int test_periodic_translation(const bool verbose){

   int ec = 0;

   const double box = 100.0;
   const std::vector<point2_t> sites = {
      point2_t(10.0,10.0), point2_t(70.0,20.0), point2_t(40.0,80.0), point2_t(85.0,85.0), point2_t(20.0,60.0)
   };
   const std::vector<double> weights(sites.size(), 0.0);

   const std::vector<polygon_t> cells_a = ::grains::internal::build_power_cells(
      sites, weights, 0.0, 0.0, box, box, true);

   std::vector<point2_t> sites_shifted(sites.size());
   for(size_t i=0; i<sites.size(); i++) sites_shifted[i] = point2_t(sites[i].x+box, sites[i].y);

   const std::vector<polygon_t> cells_b = ::grains::internal::build_power_cells(
      sites_shifted, weights, box, 0.0, 2.0*box, box, true);

   for(size_t i=0; i<sites.size(); i++){

      if(cells_a[i].size() != cells_b[i].size()){
         if(verbose) std::cout << "FAIL: site " << i << " vertex count changed under translation ("
                                << cells_a[i].size() << " vs " << cells_b[i].size() << ")" << std::endl;
         ec++;
         continue;
      }

      for(size_t v=0; v<cells_a[i].size(); v++){
         if(!nearly(cells_b[i][v].x, cells_a[i][v].x+box, 1.0e-8) || !nearly(cells_b[i][v].y, cells_a[i][v].y, 1.0e-8)){
            if(verbose) std::cout << "FAIL: site " << i << " vertex " << v << " = ("
                                   << cells_b[i][v].x << "," << cells_b[i][v].y << "), expected ("
                                   << cells_a[i][v].x+box << "," << cells_a[i][v].y << ")" << std::endl;
            ec++;
         }
      }
   }

   return ec;

}

//------------------------------------------------------------------------------
// Function to test create module power cell construction
//------------------------------------------------------------------------------
int test_grain_tessellation(const bool verbose){

   if(verbose) std::cout << "Testing grains::internal::build_power_cells()" << std::endl;

   int ec = 0;

   ec += test_square_lattice(verbose);
   ec += test_touching_discs(verbose);
   ec += test_random_configuration(verbose);
   ec += test_periodic_partition_of_unity(verbose);
   ec += test_periodic_translation(verbose);

   return ec;

}

}
}
