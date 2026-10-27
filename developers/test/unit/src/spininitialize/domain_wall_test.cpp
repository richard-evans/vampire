//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) R F L Evans 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
//   Unit tests for the domain wall spin texture, including user-defined end
//   states and neel/bloch rotation (src/spininitialize/domain_wall.cpp)
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <string>

// include header for test functions
#include "spininitialize/internal.hpp"
#include "domain_wall_test.hpp"

namespace ut{
   namespace spininitialize{

namespace si = ::spininitialize::internal;

namespace{

   const double tolerance = 1.0e-9;

   int vecerror(const double v[3], const double expected[3], const std::string test){
      for(int n = 0; n < 3; ++n){
         if(std::fabs(v[n] - expected[n]) > tolerance){
            std::cout << "FAIL: " << test << ": (" << v[0] << ", " << v[1] << ", " << v[2]
                      << ") should be (" << expected[0] << ", " << expected[1] << ", " << expected[2] << ")" << std::endl;
            return 1;
         }
      }
      return 0;
   }

   // domain wall material along axis, centred at 0.5 with a narrow width so
   // the system edges are deep in the end domains
   si::mp_t wall(const int axis){
      si::mp_t mat;
      mat.texture = si::domain_wall;
      mat.axis = axis;
      mat.centre[0] = 0.5;
      mat.width = 0.01;
      return mat;
   }

   void set_ends(si::mp_t& mat, const double s0, const double s1, const double s2, const double e0, const double e1, const double e2){
      const double s[3] = {s0, s1, s2};
      const double e[3] = {e0, e1, e2};
      for(int i = 0; i < 3; i++){
         mat.dw_start[i] = s[i];
         mat.dw_end[i] = e[i];
      }
      mat.dw_magnetisation_set = true;
   }

   void set_type(si::mp_t& mat, const si::domain_wall_type_t type, const int rotation){
      mat.dw_type = type;
      mat.dw_rotation = rotation;
      mat.dw_type_set = true;
   }

   // spin at fraction f along the wall axis
   void spin_at(const si::mp_t& mat, const double f, double s[3]){
      double pos[3] = {0.3, 0.3, 0.3};
      pos[mat.axis] = f;
      si::domain_wall_spin(mat, pos[0], pos[1], pos[2], s[0], s[1], s[2]);
   }

   // check start (f=0), centre (f=0.5) and end (f=1) spins of a wall
   int check_wall(si::mp_t mat, const bool expect_180, const double start[3], const double centre[3], const double end[3], const std::string test){

      int ec = 0;

      if(si::set_domain_wall_rotation(mat) != expect_180){
         std::cout << "FAIL: " << test << ": wrong 180 degree classification" << std::endl;
         ec++;
      }

      double s[3];
      spin_at(mat, 0.0, s); ec += vecerror(s, start,  test + " start");
      spin_at(mat, 0.5, s); ec += vecerror(s, centre, test + " centre");
      spin_at(mat, 1.0, s); ec += vecerror(s, end,    test + " end");

      return ec;

   }

}

int test_domain_wall(const bool verbose){

   int ec = 0;

   const double r = 1.0 / std::sqrt(2.0);

   const double px[3] = { 1, 0, 0}, mx[3] = {-1, 0, 0};
   const double py[3] = { 0, 1, 0}, my[3] = { 0,-1, 0};
   const double pz[3] = { 0, 0, 1}, mz[3] = { 0, 0,-1};

   //---------------------------------------------------------------------------
   // default walls: +axis -> -axis through +(axis+2)%3, as before end states
   // could be set
   //---------------------------------------------------------------------------
   ec += check_wall(wall(0), true, px, pz, mx, "default x");
   ec += check_wall(wall(1), true, py, px, my, "default y");
   ec += check_wall(wall(2), true, pz, py, mz, "default z");

   // the default wall must be unchanged by a resolve with default settings
   {
      si::mp_t mat = wall(0);
      const si::mp_t unresolved = mat;
      si::set_domain_wall_rotation(mat);
      double a[3], b[3];
      spin_at(unresolved, 0.47, a);
      spin_at(mat, 0.47, b);
      ec += vecerror(b, a, "default x matches unresolved profile");
   }

   //---------------------------------------------------------------------------
   // -z -> +z wall along x: neel and bloch with both rotation senses
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(0);
      set_ends(mat, 0, 0, -1, 0, 0, 1);
      ec += check_wall(mat, true, mz, px, pz, "-z to +z along x, default (neel +1)");
      set_type(mat, si::domain_wall_type_t::neel, -1);
      ec += check_wall(mat, true, mz, mx, pz, "-z to +z along x, neel -1");
      set_type(mat, si::domain_wall_type_t::bloch, 1);
      ec += check_wall(mat, true, mz, py, pz, "-z to +z along x, bloch +1");
      set_type(mat, si::domain_wall_type_t::bloch, -1);
      ec += check_wall(mat, true, mz, my, pz, "-z to +z along x, bloch -1");
   }

   //---------------------------------------------------------------------------
   // +y -> -y wall along x
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(0);
      set_ends(mat, 0, 1, 0, 0, -1, 0);
      ec += check_wall(mat, true, py, px, my, "+y to -y along x, neel");
      set_type(mat, si::domain_wall_type_t::bloch, 1);
      ec += check_wall(mat, true, py, pz, my, "+y to -y along x, bloch");
   }

   //---------------------------------------------------------------------------
   // tilted end states: neel centre is the axis with its m component removed
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(0);
      set_ends(mat, r, 0, r, -r, 0, -r);
      const double s[3] = { r, 0, r};
      const double c[3] = { r, 0,-r};
      const double e[3] = {-r, 0,-r};
      ec += check_wall(mat, true, s, c, e, "tilted x-z along x, neel");
   }

   //---------------------------------------------------------------------------
   // end states along the wall axis: bloch is undefined, falls back to the
   // default rotation through (axis+2)%3
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(0);
      set_ends(mat, 1, 0, 0, -1, 0, 0);
      set_type(mat, si::domain_wall_type_t::bloch, 1);
      ec += check_wall(mat, true, px, pz, mx, "head-to-head along x, bloch fallback");
      set_type(mat, si::domain_wall_type_t::bloch, -1);
      ec += check_wall(mat, true, px, mz, mx, "head-to-head along x, bloch -1 fallback");
   }

   //---------------------------------------------------------------------------
   // non-collinear walls rotate along the shortest arc, ignoring the type
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(0);
      set_ends(mat, 0, 0, 1, 1, 0, 0);
      const double c[3] = {r, 0, r};
      ec += check_wall(mat, false, pz, c, px, "90 degree +z to +x along x");
      set_type(mat, si::domain_wall_type_t::bloch, -1);
      ec += check_wall(mat, false, pz, c, px, "90 degree +z to +x along x, type ignored");
   }
   {
      // 120 degree wall in the x-y plane along y: centre at 60 degrees
      si::mp_t mat = wall(1);
      const double e[3] = {-0.5, std::sqrt(3.0)/2.0, 0};
      set_ends(mat, 1, 0, 0, e[0], e[1], e[2]);
      const double c[3] = {0.5, std::sqrt(3.0)/2.0, 0};
      ec += check_wall(mat, false, px, c, e, "120 degree in x-y along y");
   }

   //---------------------------------------------------------------------------
   // spins remain unit vectors through the wall
   //---------------------------------------------------------------------------
   {
      si::mp_t mat = wall(2);
      set_ends(mat, 0, 0.6, 0.8, 0.8, 0, -0.6);
      si::set_domain_wall_rotation(mat);
      for(int i = 0; i <= 100; i++){
         double s[3];
         spin_at(mat, 0.4 + 0.002*i, s);
         const double mod = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
         if(std::fabs(mod - 1.0) > tolerance){
            std::cout << "FAIL: unit length at f = " << 0.4 + 0.002*i << ": |s| = " << mod << std::endl;
            ec++;
            break;
         }
      }
   }

   if(verbose) std::cout << "   domain wall          : " << (ec == 0 ? "PASS" : "FAIL") << std::endl;

   return ec;

}

   }
}
