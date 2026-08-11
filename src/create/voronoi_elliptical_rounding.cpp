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
#include "create.hpp"

// create module headers
#include "internal.hpp"

namespace create{
namespace internal{

//------------------------------------------------------------------------------
//
//  Function to compute the in-plane scale factor applied to a voronoi grain
//  cross-section at height z, giving grains an ellipsoidal profile in z while
//  preserving the granular structure in the plane.
//
//  The profile blends linearly between a vertical grain wall and a true
//  ellipsoid according to the rounding parameter lambda:
//
//     R(z) = 1 - lambda * ( 1 - E(z) )
//
//  where E(z) is the normalised ellipsoidal profile
//
//     E(z) = sqrt( 1 - ( (z-z0) / c )^2 )
//
//  and c is the vertical semi-axis, taken separately above and below the
//  origin z0 so that the grain always tapers to zero cross-section at both
//  film surfaces for any choice of z0:
//
//     c = z0       for z <  z0   (lower half)
//     c = h - z0   for z >= z0   (upper half)
//
//  With z0 = h/2 the two semi-axes are equal and the profile is a true
//  ellipsoid. Limits:
//
//     lambda = 0  ->  R(z) = 1 everywhere   (vertical grain walls, as before)
//     lambda = 1  ->  R(z) = E(z)           (fully ellipsoidal grains)
//
//  At z = z0 the factor is always exactly 1, so the grains tile the plane at
//  the origin height with no volume compensation applied.
//
//------------------------------------------------------------------------------
double elliptical_rounding_factor(const double z){

   const double lambda = create::internal::voronoi_elliptical_rounding;

   // no rounding requested - return unmodified grain cross-section
   if(lambda <= 0.0) return 1.0;

   const double ssz = cs::system_dimensions[2];
   const double z0  = create::internal::voronoi_elliptical_rounding_height * ssz;

   // select the semi-axis for the relevant half of the grain
   const double c = (z >= z0) ? (ssz - z0) : z0;

   // guard against a degenerate half when the origin sits on a surface
   if(c < 1.0e-9) return 1.0 - lambda;

   const double dz  = (z - z0)/c;
   const double dz2 = dz*dz;

   // atoms beyond the poles of the ellipsoid (only reachable through rounding
   // errors at the surfaces, or for atoms slightly outside the nominal system
   // height) sit at zero ellipsoidal radius
   if(dz2 >= 1.0) return 1.0 - lambda;

   return 1.0 - lambda*( 1.0 - sqrt(1.0 - dz2) );

}

} // end of namespace internal
} // end of namespace create
