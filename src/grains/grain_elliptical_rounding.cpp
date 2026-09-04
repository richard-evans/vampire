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
#include <cmath>

// Vampire headers
#include "grains.hpp"

// grains module headers
#include "internal.hpp"

namespace grains{
namespace internal{

//------------------------------------------------------------------------------
// In-plane scale factor for a grain cross-section at height z, giving grains
// an ellipsoidal profile in z. Blends linearly between a vertical wall and a
// true ellipsoid by the rounding parameter lambda:
//    R(z) = 1 - lambda * ( 1 - E(z) ),   E(z) = sqrt( 1 - ((z-z0)/c)^2 )
// with c the vertical semi-axis above/below the origin z0 (so the grain
// tapers to zero at both film surfaces). lambda=0 gives straight walls,
// lambda=1 a full ellipsoid.
//------------------------------------------------------------------------------
double elliptical_rounding_factor(const double z){

   const double lambda = voronoi_elliptical_rounding;

   // no rounding requested - return unmodified grain cross-section
   if(lambda <= 0.0) return 1.0;

   const double ssz = cs::system_dimensions[2];
   const double z0  = voronoi_elliptical_rounding_height * ssz;

   // select the semi-axis for the relevant half of the grain
   const double c = (z >= z0) ? (ssz - z0) : z0;

   // guard against a degenerate half when the origin sits on a surface
   if(c < 1.0e-9) return 1.0 - lambda;

   const double dz  = (z - z0)/c;
   const double dz2 = dz*dz;

   // beyond the ellipsoid's poles: zero ellipsoidal radius
   if(dz2 >= 1.0) return 1.0 - lambda;

   return 1.0 - lambda*( 1.0 - sqrt(1.0 - dz2) );

}

} // end of namespace internal
} // end of namespace grains
