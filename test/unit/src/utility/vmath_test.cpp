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

// C++ standard library headers
#include <iostream>
#include <cmath>

// include header for test functions
#include "vmath.hpp"
#include "vmath_test.hpp"

namespace ut{

   // shared floating point comparison helper defined in units_test.cpp
   int floaterror(const double value, const double expected_value, const double precision, const std::string function);

   namespace utility{

int test_vmath(const bool verbose){

   if(verbose) std::cout << "Testing vmath module" << std::endl;

   int ec = 0;

   //-----------------------------------------------------------------------------------------
   // diamond of half-width 10, centred on the origin
   //-----------------------------------------------------------------------------------------
   double polyX[4] = { 10.0,  0.0, -10.0,   0.0 };
   double polyY[4] = {  0.0, 10.0,   0.0, -10.0 };

   //-----------------------------------------------------------------------------------------
   // unit factor must reproduce plain point_in_polygon() for a set of test points
   //-----------------------------------------------------------------------------------------
   const double test_x[] = { 0.0, 5.0, -5.0, 3.0, 9.9, 11.0, -9.9, 0.0 };
   const double test_y[] = { 0.0, 0.0,  0.0, 4.0, 0.0,  0.0,  0.0, 9.9 };
   for(int i=0; i<8; i++){
      const bool expected = vmath::point_in_polygon(test_x[i], test_y[i], polyX, polyY, 4);
      const bool actual   = vmath::point_in_polygon_scaled(test_x[i], test_y[i], 1.0, polyX, polyY, 4);
      if(expected != actual){
         std::cout << "FAIL: point_in_polygon_scaled(factor=1.0) disagrees with point_in_polygon at ("
                    << test_x[i] << "," << test_y[i] << ")" << std::endl;
         ec++;
      }
   }

   //-----------------------------------------------------------------------------------------
   // linear (not quadratic) scaling: half-width must scale exactly as f, not f^2
   //-----------------------------------------------------------------------------------------
   const double factors[] = { 1.0, 0.75, 0.5, 0.25 };
   const double eps = 1.0e-6;
   for(int i=0; i<4; i++){
      const double f = factors[i];
      const double expected_half_width = 10.0*f;

      // just inside the scaled diamond along +x
      const bool inside  = vmath::point_in_polygon_scaled(expected_half_width - eps, 0.0, f, polyX, polyY, 4);
      // just outside the scaled diamond along +x
      const bool outside = vmath::point_in_polygon_scaled(expected_half_width + eps, 0.0, f, polyX, polyY, 4);

      if(!inside){
         std::cout << "FAIL: point_in_polygon_scaled(factor=" << f << ") excludes point just inside the expected boundary" << std::endl;
         ec++;
      }
      if(outside){
         std::cout << "FAIL: point_in_polygon_scaled(factor=" << f << ") includes point just outside the expected boundary" << std::endl;
         ec++;
      }
   }

   //-----------------------------------------------------------------------------------------
   // degenerate and inverted scale factors contain nothing
   //-----------------------------------------------------------------------------------------
   if(vmath::point_in_polygon_scaled(0.0, 0.0, 0.0, polyX, polyY, 4) != false){
      std::cout << "FAIL: point_in_polygon_scaled(factor=0.0) should contain no points" << std::endl;
      ec++;
   }
   if(vmath::point_in_polygon_scaled(0.0, 0.0, -1.0, polyX, polyY, 4) != false){
      std::cout << "FAIL: point_in_polygon_scaled(factor<0.0) should contain no points" << std::endl;
      ec++;
   }

   //-----------------------------------------------------------------------------------------
   // repeated evaluation of a point on the scaled boundary must be deterministic
   //-----------------------------------------------------------------------------------------
   const bool boundary_first  = vmath::point_in_polygon_scaled(5.0, 5.0, 0.5, polyX, polyY, 4);
   const bool boundary_second = vmath::point_in_polygon_scaled(5.0, 5.0, 0.5, polyX, polyY, 4);
   if(boundary_first != boundary_second){
      std::cout << "FAIL: point_in_polygon_scaled is non-deterministic for a point on the scaled boundary" << std::endl;
      ec++;
   }

   if(verbose && ec == 0) std::cout << "   vmath tests passed" << std::endl;

   return ec;

}

}
}
