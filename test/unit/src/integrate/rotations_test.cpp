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
//   Unit tests for the rotation helpers used by the llg-rkmk2, llg-rkmk4 and
//   llg-dm integrators (src/integrate/rotations.hpp)
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

// include header for test functions
#include "integrate/rotations.hpp"
#include "rotations_test.hpp"

namespace ut{
   namespace integrate{

namespace rot = ::integrate::internal::rotations;

namespace{

   // check |value - expected| < tolerance and report failures
   int abserror(const double value, const double expected, const double tolerance, const std::string test){
      if(std::fabs(value - expected) < tolerance) return 0;
      std::cout << "FAIL: " << test << ": value " << value << " should be " << expected << " to within " << tolerance << std::endl;
      return 1;
   }

   int vecerror(const double v[3], const double expected[3], const double tolerance, const std::string test){
      int ec = 0;
      for(int n = 0; n < 3; ++n) ec += abserror(v[n], expected[n], tolerance, test);
      return ec > 0 ? 1 : 0;
   }

   // short scientific format for test names
   std::string sci(const double x){
      std::ostringstream ss;
      ss << x;
      return ss.str();
   }

   // random vector with components uniform in [-scale, scale]
   void random_vector(std::mt19937& gen, const double scale, double v[3]){
      std::uniform_real_distribution<double> dist(-scale, scale);
      for(int n = 0; n < 3; ++n) v[n] = dist(gen);
   }

   void random_unit_vector(std::mt19937& gen, double v[3]){
      random_vector(gen, 1.0, v);
      const double norm = std::sqrt(rot::norm_squared3(v));
      for(int n = 0; n < 3; ++n) v[n] /= norm;
   }

   // direction phi/|phi| scaled to angle th
   void scale_to_angle(const double dir[3], const double th, double phi[3]){
      const double norm = std::sqrt(rot::norm_squared3(dir));
      for(int n = 0; n < 3; ++n) phi[n] = th * dir[n] / norm;
   }

   // closed form Rodrigues rotation in long double (reference for series branch)
   void rodrigues_reference(const double phi[3], const double s[3], double out[3]){
      const long double p[3] = {phi[0], phi[1], phi[2]};
      const long double v[3] = {s[0], s[1], s[2]};
      const long double th2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
      const long double th = std::sqrt(th2);
      const long double a = std::sin(th) / th;
      const long double b = (1.0L - std::cos(th)) / th2;
      const long double c1[3] = { p[1]*v[2] - p[2]*v[1], p[2]*v[0] - p[0]*v[2], p[0]*v[1] - p[1]*v[0] };
      const long double c2[3] = { p[1]*c1[2] - p[2]*c1[1], p[2]*c1[0] - p[0]*c1[2], p[0]*c1[1] - p[1]*c1[0] };
      for(int n = 0; n < 3; ++n) out[n] = static_cast<double>(v[n] + a*c1[n] + b*c2[n]);
   }

   // closed form dexp^-1 in long double (reference for series branch)
   void dexp_inv_reference(const double phi[3], const double v_in[3], double out[3]){
      const long double p[3] = {phi[0], phi[1], phi[2]};
      const long double v[3] = {v_in[0], v_in[1], v_in[2]};
      const long double th2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
      const long double half = 0.5L * std::sqrt(th2);
      const long double beta = (1.0L - half * std::cos(half) / std::sin(half)) / th2;
      const long double c1[3] = { p[1]*v[2] - p[2]*v[1], p[2]*v[0] - p[0]*v[2], p[0]*v[1] - p[1]*v[0] };
      const long double c2[3] = { p[1]*c1[2] - p[2]*c1[1], p[2]*c1[0] - p[0]*c1[2], p[0]*c1[1] - p[1]*c1[0] };
      for(int n = 0; n < 3; ++n) out[n] = static_cast<double>(v[n] - 0.5L*c1[n] + beta*c2[n]);
   }

   // derivative of the exponential map on so(3)
   //    dexp(phi) v = v + ((1 - cos th)/th^2) phi x v + ((th - sin th)/th^3) phi x (phi x v)
   void dexp_so3(const double phi[3], const double v[3], double out[3]){
      const double th2 = rot::norm_squared3(phi);
      const double th = std::sqrt(th2);
      const double a = (1.0 - std::cos(th)) / th2;
      const double b = (th - std::sin(th)) / (th2 * th);
      double c1[3];
      rot::cross_product3(phi, v, c1);
      double c2[3];
      rot::cross_product3(phi, c1, c2);
      for(int n = 0; n < 3; ++n) out[n] = v[n] + a*c1[n] + b*c2[n];
   }

}

int test_rotations(const bool verbose){

   if(verbose) std::cout << "Testing rotation helpers for geometric LLG integrators" << std::endl;

   int ec = 0;

   std::mt19937 gen(12345);

   //---------------------------------------------------------------------------
   // Rodrigues rotation about z by pi/2 takes x to y
   //---------------------------------------------------------------------------
   {
      const double phi[3] = {0.0, 0.0, 0.5 * M_PI};
      const double s[3] = {1.0, 0.0, 0.0};
      const double expected[3] = {0.0, 1.0, 0.0};
      double out[3];
      rot::rodrigues_rotate(phi, s, out);
      ec += vecerror(out, expected, 1.0e-15, "rodrigues_rotate pi/2 about z");
   }

   //---------------------------------------------------------------------------
   // Rodrigues rotation preserves the spin length for small and large angles
   //---------------------------------------------------------------------------
   {
      const double angles[] = {1.0e-12, 1.0e-6, 0.99e-4, 1.01e-4, 1.0e-2, 0.5, 2.0, 3.0, 10.0};
      for(const double th : angles){
         for(int i = 0; i < 100; ++i){
            double dir[3], phi[3], s[3], out[3];
            random_vector(gen, 1.0, dir);
            scale_to_angle(dir, th, phi);
            random_unit_vector(gen, s);
            rot::rodrigues_rotate(phi, s, out);
            ec += abserror(std::sqrt(rot::norm_squared3(out)), 1.0, 1.0e-15, "rodrigues_rotate preserves |S| at angle " + sci(th));
            // rotation leaves the component along the rotation axis unchanged
            double axis[3];
            scale_to_angle(dir, 1.0, axis);
            ec += abserror(rot::dot3(out, axis), rot::dot3(s, axis), 1.0e-15, "rodrigues_rotate preserves S along axis at angle " + sci(th));
         }
      }
   }

   //---------------------------------------------------------------------------
   // Series (th^2 < 1e-8) and closed form branches agree near the switch
   //---------------------------------------------------------------------------
   {
      const double th2_values[] = {0.5e-8, 0.99e-8, 0.999999e-8, 1.000001e-8, 1.01e-8, 2.0e-8};
      for(const double th2 : th2_values){
         for(int i = 0; i < 100; ++i){
            double dir[3], phi[3], s[3], out[3], ref[3];
            random_vector(gen, 1.0, dir);
            scale_to_angle(dir, std::sqrt(th2), phi);
            random_unit_vector(gen, s);

            rot::rodrigues_rotate(phi, s, out);
            rodrigues_reference(phi, s, ref);
            ec += vecerror(out, ref, 1.0e-15, "rodrigues_rotate series/closed form at th^2 = " + sci(th2));

            rot::dexp_inv_so3(phi, s, out);
            dexp_inv_reference(phi, s, ref);
            ec += vecerror(out, ref, 1.0e-15, "dexp_inv_so3 series/closed form at th^2 = " + sci(th2));
         }
      }
   }

   //---------------------------------------------------------------------------
   // dexp_inv_so3 inverts dexp_so3 to round-off
   //---------------------------------------------------------------------------
   {
      const double angles[] = {1.0e-3, 1.0e-2, 0.1, 0.5, 1.0, 2.0, 3.0};
      for(const double th : angles){
         for(int i = 0; i < 100; ++i){
            double dir[3], phi[3], v[3], w[3], out[3];
            random_vector(gen, 1.0, dir);
            scale_to_angle(dir, th, phi);
            random_vector(gen, 1.0, v);
            dexp_so3(phi, v, w);
            rot::dexp_inv_so3(phi, w, out);
            ec += vecerror(out, v, 1.0e-13, "dexp_inv_so3(dexp_so3(v)) = v at angle " + sci(th));
         }
      }
   }

   //---------------------------------------------------------------------------
   // omega x S reproduces the LLG derivative used in the Heun integrator
   // (src/integrate/llg_heun.cpp) with one_oneplusalpha_sq = -gamma/(1+alpha^2)
   //---------------------------------------------------------------------------
   {
      std::uniform_real_distribution<double> alpha_dist(0.0, 1.0);
      for(int i = 0; i < 1000; ++i){
         double S[3], H[3];
         random_unit_vector(gen, S);
         random_vector(gen, 10.0, H);
         const double alpha = alpha_dist(gen);

         const double one_oneplusalpha_sq = -1.0 / (1.0 + alpha * alpha);
         const double alpha_oneplusalpha_sq = alpha * one_oneplusalpha_sq;

         double xyz[3];
         xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
         xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
         xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

         double omega[3], dS[3];
         rot::omega_llg(S, H, -one_oneplusalpha_sq, alpha, omega);
         rot::cross_product3(omega, S, dS);

         ec += vecerror(dS, xyz, 1.0e-13, "omega_llg x S equals Heun LLG derivative");
      }
   }

   if(verbose){
      if(ec == 0) std::cout << "   rotations            : PASS" << std::endl;
      else        std::cout << "   rotations            : FAIL " << ec << std::endl;
   }

   return ec;

}

} // end of integrate namespace
} // end of ut namespace
