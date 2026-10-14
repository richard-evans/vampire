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
// a flat columnar body with a rounded cap only at the top - the shape of a
// sputtered/columnar grain, rather than a lens tapering at both the substrate
// and the free surface.
//
// Grain thickness is measured against grain_film_height (falling back to the
// full simulation cell height cs::system_dimensions[2] if unset), NOT the raw
// atom z-coordinate range, so that a non-magnetic layer stacked above the
// grains (e.g. a dense dipole-field sensor slab) is left untouched.
//
//    z < z0            : R(z) = 1                                  (flat column)
//    z0 <= z < film_top : R(z) = 1 - lambda*(1 - E(z)), E(z) = sqrt(1-((z-z0)/c)^2)
//    z >= film_top      : R(z) = 1                                  (above the grain - untouched)
//
// with z0 = voronoi_elliptical_rounding_height * film_top the cap start
// height and c = film_top - z0 the cap's vertical extent. lambda=0 leaves the
// cap flat (straight/"shear" sides, i.e. no rounding at all); lambda=1 gives a
// fully rounded/ellipsoidal dome spanning the cap region.
//------------------------------------------------------------------------------
double elliptical_rounding_factor(const double z){

   const double lambda = voronoi_elliptical_rounding;

   // no rounding requested - return unmodified grain cross-section
   if(lambda <= 0.0) return 1.0;

   const double film_top = (grain_film_height > 0.0) ? grain_film_height : cs::system_dimensions[2];
   const double z0 = voronoi_elliptical_rounding_height * film_top;

   // below the cap start: flat column, full width
   if(z < z0) return 1.0;

   // above the grain film (e.g. inside a sensor/fill layer stacked on top): untouched
   if(z >= film_top) return 1.0;

   const double c = film_top - z0;

   // guard against a degenerate cap when the start height sits on the film top
   if(c < 1.0e-9) return 1.0 - lambda;

   const double dz  = (z - z0)/c;
   const double dz2 = dz*dz;

   // beyond the dome's pole: zero cap radius
   if(dz2 >= 1.0) return 1.0 - lambda;

   return 1.0 - lambda*( 1.0 - sqrt(1.0 - dz2) );

}

} // end of namespace internal
} // end of namespace grains
