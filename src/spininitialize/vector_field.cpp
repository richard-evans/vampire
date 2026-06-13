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

         const std::vector<field_point_t>& points = vector_field_data[mat.vector_field_id];

         const double power = 2.0; // inverse-distance weighting exponent
         const double tol   = 1.0e-12;

         double sum_w  = 0.0;
         double sum_mx = 0.0;
         double sum_my = 0.0;
         double sum_mz = 0.0;

         for(size_t i = 0; i < points.size(); i++){

            const double dx = fx - points[i].x;
            const double dy = fy - points[i].y;
            const double dz = fz - points[i].z;

            const double d2 = dx*dx + dy*dy + dz*dz;

            // if the atom lies (almost) exactly on a field point, use it directly
            if(d2 < tol){
               sx = points[i].mx;
               sy = points[i].my;
               sz = points[i].mz;
               return;
            }

            const double w = 1.0 / pow(d2, 0.5*power);

            sum_w  += w;
            sum_mx += w * points[i].mx;
            sum_my += w * points[i].my;
            sum_mz += w * points[i].mz;

         }

         sx = sum_mx / sum_w;
         sy = sum_my / sum_w;
         sz = sum_mz / sum_w;

      }

   } // end of internal namespace

} // end of spininitialize namespace
