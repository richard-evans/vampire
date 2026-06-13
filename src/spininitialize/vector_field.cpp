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

// Vampire headers
#include "spininitialize.hpp"

// spininitialize module headers
#include "internal.hpp"

namespace spininitialize{

   namespace internal{

      //-----------------------------------------------------------------------------
      // vector_field: spin direction interpolated from a user-supplied list of
      // (x,y,z,mx,my,mz) points (fractional coordinates and unit vectors) using
      // inverse-distance weighting over all points in the field.
      //
      // If the atom coincides (to within a small tolerance) with one of the
      // supplied points, that point's direction is used directly.
      //-----------------------------------------------------------------------------
      void vector_field_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         // look up the list of (x,y,z,mx,my,mz) points loaded from the field file
         // for this material (cached in vector_field_data, see parse.cpp)
         const std::vector<field_point_t>& points = vector_field_data[mat.vector_field_id];

         // inverse-distance weighting (IDW): each field point contributes a
         // weight w_i = 1 / d_i^power, where d_i is the distance from the atom
         // to that point. power = 2 gives an "inverse-square" weighting, so
         // nearby points dominate strongly over distant ones.
         const double power = 2.0; // inverse-distance weighting exponent
         const double tol   = 1.0e-12; // distance^2 below which a point is treated as coincident

         // running totals for the weighted sum: sum_w accumulates the
         // normalisation (sum of all weights), sum_m{x,y,z} accumulate the
         // weighted vector components
         double sum_w  = 0.0;
         double sum_mx = 0.0;
         double sum_my = 0.0;
         double sum_mz = 0.0;

         // loop over every point in the field and accumulate its contribution
         for(size_t i = 0; i < points.size(); i++){

            // displacement (in fractional coordinates) from the field point to this atom
            const double dx = fx - points[i].x;
            const double dy = fy - points[i].y;
            const double dz = fz - points[i].z;

            // squared distance (avoids an unnecessary sqrt in the common case)
            const double d2 = dx*dx + dy*dy + dz*dz;

            // if the atom lies (almost) exactly on a field point, use its
            // direction directly to avoid dividing by (near) zero
            if(d2 < tol){
               sx = points[i].mx;
               sy = points[i].my;
               sz = points[i].mz;
               return;
            }

            // weight w_i = d_i^(-power) = (d2)^(-power/2)
            const double w = 1.0 / pow(d2, 0.5*power);

            sum_w  += w;
            sum_mx += w * points[i].mx;
            sum_my += w * points[i].my;
            sum_mz += w * points[i].mz;

         }

         // weighted average direction: sum(w_i * m_i) / sum(w_i)
         // this is not generally a unit vector, so it is renormalised by the
         // caller (see get_spin_direction in textures.cpp)
         sx = sum_mx / sum_w;
         sy = sum_my / sum_w;
         sz = sum_mz / sum_w;

      }

   } // end of internal namespace

} // end of spininitialize namespace
