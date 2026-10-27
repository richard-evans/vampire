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

      namespace {

         const double pi = 3.14159265358979323846;

         // end states closer than this to antiparallel (1 + m0.m1 < tol)
         // are treated as a 180 degree wall
         const double antiparallel_tolerance = 1.0e-10;

         // normalise a 3-vector in place; returns false for a zero vector
         bool unit(double v[3]){
            const double mod = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            if(mod < 1.0e-12) return false;
            for(int i = 0; i < 3; i++) v[i] /= mod;
            return true;
         }

      }

      //-----------------------------------------------------------------------------
      // set_domain_wall_rotation: resolve the rotation of a domain wall into
      // an orthonormal pair (dw_m0, dw_m1) and a total rotation angle
      // dw_angle, so that domain_wall_spin() only needs to evaluate
      //    m = cos(phi) dw_m0 + sin(phi) dw_m1
      //
      // dw_m0 is the start magnetisation. For a 180 degree wall dw_m1 is the
      // spin direction at the wall centre, chosen by dw_type and dw_rotation
      // (see domain_wall_type_t). If that direction is undefined (end states
      // along the wall axis) the wall rotates through +/- the (axis+2)%3
      // direction, which reproduces the original +axis -> -axis wall. For a
      // non-collinear wall dw_m1 is the component of the end state
      // perpendicular to the start, and the type is ignored.
      //
      // Returns true if the wall is a 180 degree wall.
      //-----------------------------------------------------------------------------
      bool set_domain_wall_rotation(mp_t& mat){

         // end states, defaulting to +axis -> -axis
         double m0[3] = {0.0, 0.0, 0.0};
         double m1[3] = {0.0, 0.0, 0.0};
         if(mat.dw_magnetisation_set){
            for(int i = 0; i < 3; i++){
               m0[i] = mat.dw_start[i];
               m1[i] = mat.dw_end[i];
            }
         }
         else{
            m0[mat.axis] =  1.0;
            m1[mat.axis] = -1.0;
         }
         unit(m0);
         unit(m1);

         for(int i = 0; i < 3; i++) mat.dw_m0[i] = m0[i];

         const double c = m0[0]*m1[0] + m0[1]*m1[1] + m0[2]*m1[2];

         //---------------------------------------------------------------
         // non-collinear wall: rotate along the shortest arc from m0 to m1
         //---------------------------------------------------------------
         if(1.0 + c > antiparallel_tolerance){
            for(int i = 0; i < 3; i++) mat.dw_m1[i] = m1[i] - c * m0[i];
            unit(mat.dw_m1);
            mat.dw_angle = std::acos(std::fmin(1.0, std::fmax(-1.0, c)));
            return false;
         }

         //---------------------------------------------------------------
         // 180 degree wall: centre direction set by the wall type
         //---------------------------------------------------------------
         double a[3] = {0.0, 0.0, 0.0}; // unit vector along the wall axis
         a[mat.axis] = 1.0;

         double centre[3];
         if(mat.dw_type == domain_wall_type_t::neel){
            // axis with its component along m0 removed
            const double am = a[0]*m0[0] + a[1]*m0[1] + a[2]*m0[2];
            for(int i = 0; i < 3; i++) centre[i] = a[i] - am * m0[i];
         }
         else{
            // a x m0
            centre[0] = a[1]*m0[2] - a[2]*m0[1];
            centre[1] = a[2]*m0[0] - a[0]*m0[2];
            centre[2] = a[0]*m0[1] - a[1]*m0[0];
         }

         // end states along the wall axis: neither neel nor bloch is
         // defined, so rotate through the (axis+2)%3 direction, which is
         // perpendicular to m0 in this case
         if(!unit(centre)){
            for(int i = 0; i < 3; i++) centre[i] = 0.0;
            centre[(mat.axis+2)%3] = 1.0;
         }

         for(int i = 0; i < 3; i++) mat.dw_m1[i] = double(mat.dw_rotation) * centre[i];
         mat.dw_angle = pi;

         return true;

      }

      //-----------------------------------------------------------------------------
      // domain_wall: 1D tanh-profile domain wall.
      //
      // Along the chosen axis the spin rotates from dw_m0 (far below the wall
      // centre) through dw_m1 to the end state (far above the wall centre),
      // using the rotation resolved by set_domain_wall_rotation(). By default
      // this is a 180 degree wall from +<axis> through +<(axis+2)%3> to
      // -<axis>, i.e. for domain-wall-x the rotation plane is x-z, for
      // domain-wall-y it is y-x and for domain-wall-z it is z-y.
      //
      // mat.centre[0] is the position of the wall centre (fraction of system
      // size along axis) and mat.width controls the sharpness of the wall
      // (fraction of system size); smaller values give a sharper wall.
      //-----------------------------------------------------------------------------
      void domain_wall_spin(const mp_t& mat, const double fx, const double fy, const double fz, double& sx, double& sy, double& sz){

         // pick out the fractional coordinate along the wall's propagation
         // axis (0=x, 1=y, 2=z), as set by domain-wall-x/y/z
         const double pos[3] = {fx, fy, fz};

         // t is a dimensionless "distance from the wall centre" measured in
         // units of the wall width: t = 0 at the centre, |t| >> 1 far away.
         // mat.centre[0] and mat.width are both fractions of the system size.
         const double t = (pos[mat.axis] - mat.centre[0]) / mat.width;

         // standard tanh-profile domain wall: the polar angle theta follows
         //   theta(t) = 2*atan(exp(t))
         // which smoothly varies from 0 (t -> -inf) through pi/2 (t = 0, the
         // wall centre) to pi (t -> +inf).
         const double theta = 2.0 * atan(exp(t));

         // scale to the total rotation angle of the wall (exactly theta for
         // a 180 degree wall)
         const double phi = theta * (mat.dw_angle / pi);

         const double c = cos(phi);
         const double s = sin(phi);

         sx = c * mat.dw_m0[0] + s * mat.dw_m1[0];
         sy = c * mat.dw_m0[1] + s * mat.dw_m1[1];
         sz = c * mat.dw_m0[2] + s * mat.dw_m1[2];

      }

   } // end of internal namespace

} // end of spininitialize namespace
