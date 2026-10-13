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
//   This file is a derivative work of JAMS (https://github.com/stonerlab/jams),
//   ported from src/jams/solvers/llg_rkmk_functions.h
//   (version3 branch, commit 234ba11).
//   Original implementation by Dr Joseph Barker, University of Leeds.
//
//   The original JAMS code is distributed under the MIT License, reproduced
//   below in full as required:
//
//   MIT License
//
//   Copyright (c) 2009-2013 Joseph Barker (University of York, UK),
//                 2014-2018 Joseph Barker (Tohoku University, Japan),
//                 2018-2024 Joseph Barker and collaborators (University of Leeds, UK)
//
//   Permission is hereby granted, free of charge, to any person obtaining a copy
//   of this software and associated documentation files (the "Software"), to deal
//   in the Software without restriction, including without limitation the rights
//   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//   copies of the Software, and to permit persons to whom the Software is
//   furnished to do so, subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in all
//   copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.
//
//------------------------------------------------------------------------------
//

#ifndef INTEGRATE_ROTATIONS_H_
#define INTEGRATE_ROTATIONS_H_

//------------------------------------------------------------------------------
// Inline vector and rotation helpers for the geometric (rotation based) LLG
// integrators: Runge-Kutta-Munthe-Kaas (RKMK2, RKMK4) and Depondt-Mertens.
//
// In these schemes the LLG equation is written as a rotation of the spin
//
//    dS/dt = omega x S,   omega = gyro (H + alpha S x H)
//
// with gyro = gamma/(1+alpha^2). Each step rotates the spin by a rotation
// vector phi = dt omega using the Rodrigues formula, so the spin length is
// preserved exactly without renormalisation.
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>

namespace integrate{
namespace internal{
namespace rotations{

//------------------------------------------------------------------------------
// Basic 3-vector operations
//------------------------------------------------------------------------------
inline double dot3(const double a[3], const double b[3]){
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm_squared3(const double v[3]){
   return dot3(v, v);
}

inline void cross_product3(const double a[3], const double b[3], double out[3]){
   out[0] = a[1]*b[2] - a[2]*b[1];
   out[1] = a[2]*b[0] - a[0]*b[2];
   out[2] = a[0]*b[1] - a[1]*b[0];
}

//------------------------------------------------------------------------------
// Angular velocity of the LLG equation, omega = gyro (H + alpha S x H),
// such that dS/dt = omega x S
//------------------------------------------------------------------------------
inline void omega_llg(const double s[3], const double h[3], const double gyro, const double alpha, double out[3]){
   out[0] = gyro * (h[0] + alpha * (s[1]*h[2] - s[2]*h[1]));
   out[1] = gyro * (h[1] + alpha * (s[2]*h[0] - s[0]*h[2]));
   out[2] = gyro * (h[2] + alpha * (s[0]*h[1] - s[1]*h[0]));
}

//------------------------------------------------------------------------------
// Rodrigues rotation of s by rotation vector phi (angle |phi| about phi/|phi|)
//
//    out = s + (sin th / th) phi x s + ((1 - cos th) / th^2) phi x (phi x s)
//
// A Taylor series is used for small angles (th^2 < 1e-8) to avoid
// cancellation errors.
//------------------------------------------------------------------------------
inline void rodrigues_rotate(const double phi[3], const double s[3], double out[3]){

   const double th2 = norm_squared3(phi);

   double c1[3];
   cross_product3(phi, s, c1);
   double c2[3];
   cross_product3(phi, c1, c2);

   if(th2 < 1e-8){
      const double th4 = th2 * th2;
      const double a = 1.0 - th2 * (1.0 / 6.0) + th4 * (1.0 / 120.0);
      const double b = 0.5 - th2 * (1.0 / 24.0) + th4 * (1.0 / 720.0);

      for(int n = 0; n < 3; ++n) out[n] = s[n] + a * c1[n] + b * c2[n];
      return;
   }

   const double th = std::sqrt(th2);
   const double sin_th = std::sin(th);
   const double cos_th = std::cos(th);

   const double a = sin_th / th;
   const double b = (1.0 - cos_th) / th2;

   for(int n = 0; n < 3; ++n) out[n] = s[n] + a * c1[n] + b * c2[n];

}

//------------------------------------------------------------------------------
// Inverse of the derivative of the exponential map on so(3)
//
//    out = v - 1/2 phi x v + beta phi x (phi x v)
//    beta = (1 - (th/2) cot(th/2)) / th^2
//
// A Taylor series is used for small angles (th^2 < 1e-8).
//------------------------------------------------------------------------------
inline void dexp_inv_so3(const double phi[3], const double v[3], double out[3]){

   const double th2 = norm_squared3(phi);

   double c1[3];
   cross_product3(phi, v, c1);
   double c2[3];
   cross_product3(phi, c1, c2);

   if(th2 < 1e-8){
      const double th4 = th2 * th2;
      const double beta = (1.0 / 12.0) + th2 * (1.0 / 720.0) + th4 * (1.0 / 30240.0);

      for(int n = 0; n < 3; ++n) out[n] = v[n] - 0.5 * c1[n] + beta * c2[n];
      return;
   }

   const double th = std::sqrt(th2);
   const double half = 0.5 * th;

   const double sin_half = std::sin(half);
   const double cos_half = std::cos(half);

   const double cot_half = cos_half / sin_half;
   const double beta = (1.0 / th2) * (1.0 - half * cot_half);

   for(int n = 0; n < 3; ++n) out[n] = v[n] - 0.5 * c1[n] + beta * c2[n];

}

//------------------------------------------------------------------------------
// Stochastic substep: rotate s by dt omega(s, h_noise) with the spin frozen
// at the start of the substep
//------------------------------------------------------------------------------
inline void noise_step_rodrigues(const double s[3], const double h_noise[3], const double gyro, const double alpha, const double dt, double out[3]){

   double omega[3];
   omega_llg(s, h_noise, gyro, alpha, omega);

   const double phi[3] = {dt * omega[0], dt * omega[1], dt * omega[2]};
   rodrigues_rotate(phi, s, out);

}

} // end of rotations namespace
} // end of internal namespace
} // end of integrate namespace

#endif //INTEGRATE_ROTATIONS_H_
